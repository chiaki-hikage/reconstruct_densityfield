#define pi 3.141592653589793238462643
#include<iostream>
#include<time.h>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"recon.h"
#include"fourier_ver3.h"

using namespace std;

void assign_density(double **d, double *len, const int& n, double ***f, int *npix);
void read_and_assign_density(double ***f, double *len, int *npix, const Parlist par);
double**initmat2(const int& n1,const int& n2);
double***initmat3(const int& n1,const int& n2,const int& n3);
double****initmat4(const int& n1,const int& n2,const int& n3, const int& n4);
void make_density(double ***f, int *npix, const int& nsub);
void make_norm(double ***f, int *npix, double *len, double ***pshot, const int& nsub);
void smooth_field(double ***f, int *npix, const double& s);
void smooth_field_sub(double ***f, int *npix, const double& s, const int& nsub);
void deriv_field(double ***f, double ****v, int *npix, double *len);
void deriv_field_sub(double ***f, double ****v, int *npix, double *len, const int& nsub);
void iso_shiftfield(double ****v, int *npix, double *len, const double& fz);
void reconstruct_data(double ***f, double *len, int *npix, char *inf, double ****v, const double& eps, const int& nsnap, const Parlist par);
void reconstruct_add(double **d, const int& np, double ****v, int *npix, const double& eps, double *len);
void reconstruct_grid(double ***f, double ****v, int *npix, const double& eps, double *len);
void calc_pk(double ***f, double *len, int *npix, const double& kmin, const double& kmax, const int& nkbin, const double& shotn, char* outf);
void calc_pk_sub(double ***f, double *len, int *npix, const double& kmin, const double& kmax, const int& nkbin, double ***shotn, char* outf, const int& nsub);
void calc_pk_cross(double ***f, double ***f2, double *len, int *npix, const double& kmin, const double& kmax, const int& nkbin, const double& shotn, const double& shotn2, char* outf);
void calc_pk_cross_sub(double ***f, double ***f2, double *len, int *npix, const double& kmin, const double& kmax, const int& nkbin, double ***shotn, double ***shotn2, char* outf, const int& nsub);
void ofwritefunc(char*string, double*x, double*y, const int& n);
void ofwritefunc_2cols(char*string, double*x, double**y, const int& n, const int& m);
void ofwritefunc_3cols(char*string, double*x, double***y, const int& n, const int& m, const int& l);
double norm(double *x, const int& n);
double gauss(const double& k);
double hrate(const Parlist par);

typedef struct
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int flag_stellarage;
  int flag_metals;
  int hashtabsize;
  char fill[84];		/* fills to 256 Bytes */
}
io_header;

int main(int argc, char** argv){

  Parlist par;
  par.defaultval();
  par.argument(argc,argv);

  //!!! k_min, k_max for output pk
  par.kmin=0.01/(2*pi/par.len*par.npix/2);
  par.kmax=0.81/(2*pi/par.len*par.npix/2);
  par.nkbin=80;

  char outf[100];
  double len[3]={par.len, par.len, par.len};
  double lensub[3]={len[0]/par.nsub, len[1]/par.nsub, len[2]/par.nsub};
  int npix[3]={par.npix,par.npix,par.npix};
  long int nparticles=pow(par.npart,3);
  double boxsize=pow(par.len,3);

  double ***d;
  double nden=(double)nparticles/boxsize;
  double pshot=(double)1/nden;
  double ***pshot_d=initmat3(par.nsub,par.nsub,par.nsub);
  double ***pshot_r=initmat3(par.nsub,par.nsub,par.nsub);
  double norm_f=nden*sqrt(boxsize);
  cout << "number of mass particles: " << nparticles << ", density: " << nden << ", z=" << par.z << endl;
  cout << "pshot: " << pshot << endl;

  double ***f=initmat3(npix[0],npix[1],npix[2]+2);
  double ***f_pre=initmat3(npix[0],npix[1],npix[2]+2);
  double ***pshot_d_pre=initmat3(par.nsub,par.nsub,par.nsub);
  double pixden_f=(double)nparticles/npix[0]/npix[1]/npix[2];

  cout << "read particle data and make density fields" << endl;
  read_and_assign_density(f,len,npix,par);

  if (par.rec_flag[0]=='y') {
    cout << "start reconstruction with Rs=" << par.rscale << endl;
    cout << "smooth density field" << endl;

    if (par.calc_crossP[0]=='y') {
      for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) for (int k=0;k<npix[2];k++) f_pre[i][j][k]=f[i][j][k];
      if (par.nsub>1) {
	make_norm(f_pre,npix,lensub,pshot_d_pre,par.nsub);
      } else {
	for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) for (int k=0;k<npix[2];k++) f_pre[i][j][k]=(f_pre[i][j][k]-pixden_f)/norm_f;
      }
    }

    if (par.nsub>1) {
      make_density(f,npix,par.nsub);
      smooth_field_sub(f,npix,par.rscale/lensub[0],par.nsub);
    } else {
      for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) for (int k=0;k<npix[2];k++) f[i][j][k]=f[i][j][k]/pixden_f-1.;
      smooth_field(f,npix,par.rscale/len[0]);
    }

    cout << "obtain displacement field via inverse ZA: s(q)" << endl;
    double ****v=initmat4(3,npix[0],npix[1],npix[2]+2);
    if (par.nsub>1) {
      deriv_field_sub(f,v,npix,lensub,par.nsub);
    } else {
      deriv_field(f,v,npix,len);
    }

    cout << "displace data by s(x) and obtain delta_d" << endl;
    for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) for (int k=0;k<npix[2]+2;k++) f[i][j][k]=0;
    reconstruct_data(f,len,npix,par.inf,v,1.,par.nsnap,par);
    if (par.nsub>1) {
      make_norm(f,npix,lensub,pshot_d,par.nsub);
    } else {
      for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) for (int k=0;k<npix[2];k++) f[i][j][k]=(f[i][j][k]-pixden_f)/norm_f;
    }

    if (par.space[0]=='z'&&par.rectype[0]=='i') iso_shiftfield(v,npix,len,par.fz);

    cout << "displace random and obtain delta_s" << endl;
    double ***g=initmat3(npix[0],npix[1],npix[2]+2);
    reconstruct_grid(g,v,npix,1.,len);
    delete [] ***v; delete [] **v; delete [] *v; delete [] v;

    if (par.nsub>1) {
      make_norm(g,npix,lensub,pshot_r,par.nsub);
    } else {
      double norm_g=(npix[0]/len[0])*(npix[1]/len[1])*(npix[2]/len[2])*sqrt(boxsize);
      double pixden_g=1.;
      for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) for (int k=0;k<npix[2];k++) g[i][j][k]=(g[i][j][k]-pixden_g)/norm_g;
    }
    
    cout << " obtain reconstruct field: delta_d - delta_s  " << endl;
    for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) for (int k=0;k<npix[2];k++) f[i][j][k]-=g[i][j][k];
    delete [] **g; delete [] *g; delete[] g;
    
  } else {
    cout << "No reconstruction" << endl;
    if (par.nsub>1) {
      make_norm(f,npix,lensub,pshot_d,par.nsub);
    } else {
      for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) for (int k=0;k<npix[2];k++) f[i][j][k]=(f[i][j][k]-pixden_f)/norm_f;
    }
  }

  cout << "compute power spectrum" << endl;
  cout << "output file: " << par.outf << endl;

  if (par.calc_crossP[0]=='y') {
    if (par.nsub>1) {
      calc_pk_cross_sub(f,f_pre,lensub,npix,par.kmin,par.kmax,par.nkbin,pshot_d,pshot_d_pre,par.outf,par.nsub);
    } else {
      calc_pk_cross(f,f_pre,len,npix,par.kmin,par.kmax,par.nkbin,pshot,pshot,par.outf);
    }
    delete [] **f_pre; delete [] *f_pre; delete[] f_pre;
  } else {
    if (par.nsub>1) {
      calc_pk_sub(f,lensub,npix,par.kmin,par.kmax,par.nkbin,pshot_d,par.outf,par.nsub);
    } else {
      calc_pk(f,len,npix,par.kmin,par.kmax,par.nkbin,pshot,par.outf);
    }
  }

  delete [] **f; delete [] *f; delete[] f;

  return 0;
}

void read_and_assign_density(double ***f, double *len, int *npix, const Parlist par) {
  ifstream ifs;
  int dummy;
  FILE *fd;
  io_header header;
  char inf[150];
  long int nptot;
  double ahz=100*hrate(par)/(1+par.z);

  nptot=0;
  for (int isnap=0;isnap<par.nsnap;isnap++) {
    sprintf(inf,"%s.%d",par.inf,isnap);
    ifs.open(inf,ios::in|ios::binary);
    if (!ifs) throw ios::failure("Failed to open input file");
    fd=fopen(inf,"r");
    fread(&dummy,sizeof(int),1,fd);
    fread(&header,sizeof(io_header),1,fd);
    fread(&dummy,sizeof(int),1,fd);
    cout << inf << " " << header.npart[1] << endl;
    fread(&dummy,sizeof(int),1,fd);

    int np=header.npart[1];
    double **d=initmat2(np,3);
    for (int i=0;i<np;i++) {
      float x[3]; 
      fread(x,sizeof(float),3,fd);
      for (int j=0;j<3;j++) d[i][j]=x[j];
    }

    nptot+=(long)np;   // (long)
    cout << nptot << endl;
    fread(&dummy,sizeof(int),1,fd);

    if (par.space[0]=='z') {
      float r_vel=0;
      fread(&dummy,sizeof(int),1,fd);
      for (int i=0;i<np;i++) {
	float v[3]; 
	fread(v,sizeof(float),3,fd);
	for (int j=0;j<3;j++) v[j]/=sqrt(1+par.z);
	for (int j=0;j<3;j++) v[j]*=par.fz_boost;
	d[i][2]=fmod(d[i][2]+v[2]/ahz+len[2],len[2]); 
	r_vel+=pow(v[2]/ahz,2);
      }
      fread(&dummy,sizeof(int),1,fd);
      r_vel=sqrt(r_vel/np);
    }

    fclose(fd);   // fclose
    ifs.close();
    assign_density(d,len,np,f,npix);
    delete [] *d; delete [] d;
  }
}

void reconstruct_data(double ***f, double *len, int *npix, char *infin, double ****v, const double& eps, const int& nsnap, const Parlist par) {
  ifstream ifs;
  int dummy;
  long int nptot;
  FILE *fd;
  io_header header;
  char inf[150];

  double ahz=100*hrate(par)/(1+par.z);

  nptot=0;
  for (int isnap=0;isnap<nsnap;isnap++) {
    sprintf(inf,"%s.%d",infin,isnap);
    ifs.open(inf,ios::in|ios::binary);
    if (!ifs) throw ios::failure("Failed to open input file");
    fd=fopen(inf,"r");
    fread(&dummy,sizeof(int),1,fd);
    fread(&header,sizeof(io_header),1,fd);
    fread(&dummy,sizeof(int),1,fd);
    cout << inf << " " << header.npart[1] << endl;
    fread(&dummy,sizeof(int),1,fd);

    int np=header.npart[1];
    double **d=initmat2(np,3);
    for (int i=0;i<np;i++) {
      float x[3]; 
      fread(x,sizeof(float),3,fd);
      for (int j=0;j<3;j++) d[i][j]=x[j];
    }
    cout << np << endl;
    nptot+=(long)np;           // (long)
    fread(&dummy,sizeof(int),1,fd);

    if (par.space[0]=='z') {
      float r_vel=0;
      fread(&dummy,sizeof(int),1,fd);
      for (int i=0;i<np;i++) {
	float v[3]; 
	fread(v,sizeof(float),3,fd);
	for (int j=0;j<3;j++) v[j]/=sqrt(1+par.z);
	for (int j=0;j<3;j++) v[j]*=par.fz_boost;
	d[i][2]=fmod(d[i][2]+v[2]/ahz+len[2],len[2]); 
	r_vel+=pow(v[2]/ahz,2);
      }
      fread(&dummy,sizeof(int),1,fd);
    }

    fclose(fd);    // fclose
    ifs.close();
    reconstruct_add(d,np,v,npix,eps,len);
    assign_density(d,len,np,f,npix);
    delete [] *d; delete [] d;
  }
}

void reconstruct_add(double **d, const int& np, double ****v, int *npix, const double& eps, double *len) {
  int ip[3],ip1[3],ipn[3];
  double fp[3],mean[3]={0,0,0},disp[3],maxdis=0.;
  for (int i=0;i<np;i++) {
    for (int j=0;j<3;j++) {
      ip[j]=(int)(floor(d[i][j]/len[j]*npix[j]-0.5));
      fp[j]=d[i][j]/len[j]*npix[j]-0.5-ip[j];
      ip[j]=(ip[j]+npix[j])%npix[j];
      ip1[j]=(ip[j]+1)%npix[j];
      ipn[j]=(int)(d[i][j]/len[j]*npix[j]);
      if (fp[j]<0||fp[j]>1) {
	cout << "fp is not in the range of [0,1] " << fp[j] << endl;
	throw ios::failure("error");
      }
    }
    for (int j=0;j<3;j++) {
      disp[j]=eps*v[j][ip[0]][ip[1]][ip[2]]*(1-fp[0])*(1-fp[1])*(1-fp[2]);
      disp[j]+=eps*v[j][ip1[0]][ip[1]][ip[2]]*fp[0]*(1-fp[1])*(1-fp[2]);
      disp[j]+=eps*v[j][ip[0]][ip1[1]][ip[2]]*(1-fp[0])*fp[1]*(1-fp[2]);
      disp[j]+=eps*v[j][ip[0]][ip[1]][ip1[2]]*(1-fp[0])*(1-fp[1])*fp[2];
      disp[j]+=eps*v[j][ip1[0]][ip1[1]][ip[2]]*fp[0]*fp[1]*(1-fp[2]);
      disp[j]+=eps*v[j][ip[0]][ip1[1]][ip1[2]]*(1-fp[0])*fp[1]*fp[2];
      disp[j]+=eps*v[j][ip1[0]][ip[1]][ip1[2]]*fp[0]*(1-fp[1])*fp[2];
      disp[j]+=eps*v[j][ip1[0]][ip1[1]][ip1[2]]*fp[0]*fp[1]*fp[2];
      d[i][j]+=disp[j];
      d[i][j]=fmod(d[i][j]+len[j],len[j]);
    }
    double disp_abs=sqrt(disp[0]*disp[0]+disp[1]*disp[1]+disp[2]*disp[2]);
    if (maxdis<disp_abs) maxdis=disp_abs;
    for (int j=0;j<3;j++) mean[j]+=disp[j]*disp[j];
  }
  for (int j=0;j<3;j++) mean[j]=sqrt(mean[j]/np);
  cout << "Mean displacement: " << mean[0] << " " << mean[1] << " " << mean[2] << " ";
  cout << sqrt(mean[0]*mean[0]+mean[1]*mean[1]+mean[2]*mean[2]) << endl;
  cout << "Max displacement: " << maxdis << "[Mpc/h]" << endl;
}

double**initmat2(const int& n1,const int& n2) {
  double**f=new double*[n1];
  f[0]=new double[n1*n2];
  for (int i=1;i<n1;i++) f[i]=f[0]+i*n2;
  for (int i=0;i<n1;i++) for (int j=0;j<n2;j++) f[i][j]=0;
  return f;
}

double***initmat3(const int& n1,const int& n2,const int& n3) {
  // pointers to first level
  double***f=new double**[n1];
  
  // pointers to second level
  f[0]=new double*[(long)n1*n2];     // (long)
  for (int i=1;i<n1;i++) f[i]=f[i-1]+n2;
  
  // pointers to third level
  f[0][0]=new double[(long)n1*n2*n3];    // (long)
  for (int i=1;i<n1;i++) f[i][0]=f[i-1][0]+(long)n2*n3;    // (long)
  for (int i=0;i<n1;i++) for(int j=1;j<n2;j++) f[i][j]=f[i][j-1]+n3;
  for (int i=0;i<n1;i++) for(int j=0;j<n2;j++) for(int k=0;k<n3;k++) f[i][j][k]=(double)0;
  return f;
}

double****initmat4(const int& n1,const int& n2,const int& n3,const int& n4) {
  // pointers to first level
  double****f=new double***[n1];
  
  // pointers to second level
  f[0]=new double**[(long)n1*n2];     // (long)
  for (int i=1;i<n1;i++) f[i]=f[i-1]+n2;
  
  // pointers to third level
  f[0][0]=new double*[(long)n1*n2*n3];   // (long)
  for (int i=1;i<n1;i++) f[i][0]=f[i-1][0]+(long)n2*n3;    // (long)
  for (int i=0;i<n1;i++) for(int j=1;j<n2;j++) f[i][j]=f[i][j-1]+n3;
  
  // pointers to third level
  f[0][0][0]=new double[(long)n1*n2*n3*n4];   // (long)
  for (int i=1;i<n1;i++) f[i][0][0]=f[i-1][0][0]+(long)n2*n3*n4; //  (long)
  for (int i=0;i<n1;i++) for(int j=1;j<n2;j++) f[i][j][0]=f[i][j-1][0]+(long)n3*n4;  // (long)
  for (int i=0;i<n1;i++) for(int j=0;j<n2;j++) for(int k=1;k<n3;k++) f[i][j][k]=f[i][j][k-1]+n4;
  
  for (int i=0;i<n1;i++) for(int j=0;j<n2;j++) for(int k=0;k<n3;k++) for (int l=0;l<n4;l++) f[i][j][k][l]=(double)0;
  return f;
}

void make_density(double ***f, int *npix, const int& nsub) {
  int npsub[3]={npix[0]/nsub,npix[1]/nsub,npix[2]/nsub};
  for (int isub=0;isub<nsub;isub++) for (int jsub=0;jsub<nsub;jsub++) for (int ksub=0;ksub<nsub;ksub++) {
    double pixden_s;
    pixden_s=0;
    for (int i=0;i<npsub[0];i++) for (int j=0;j<npsub[1];j++) for (int k=0;k<npsub[2];k++) {
      int ip=i+isub*npsub[0], jp=j+jsub*npsub[1], kp=k+ksub*npsub[2];
      pixden_s+=f[ip][jp][kp];
    }
    pixden_s/=npsub[0]*npsub[1]*npsub[2];
    for (int i=0;i<npsub[0];i++) for (int j=0;j<npsub[1];j++) for (int k=0;k<npsub[2];k++) {
      int ip=i+isub*npsub[0], jp=j+jsub*npsub[1], kp=k+ksub*npsub[2];
      f[ip][jp][kp]=f[ip][jp][kp]/pixden_s-1.;
    }
  }
}

void make_norm(double ***f, int *npix, double *lensub, double ***pshot, const int& nsub) {
  int npsub[3]={npix[0]/nsub,npix[1]/nsub,npix[2]/nsub};
  double subboxsize=lensub[0]*lensub[1]*lensub[2];
  for (int isub=0;isub<nsub;isub++) for (int jsub=0;jsub<nsub;jsub++) for (int ksub=0;ksub<nsub;ksub++) {
    double npartsub;
    npartsub=0.;
    for (int i=0;i<npsub[0];i++) for (int j=0;j<npsub[1];j++) for (int k=0;k<npsub[2];k++) {
      int ip=i+isub*npsub[0], jp=j+jsub*npsub[1], kp=k+ksub*npsub[2];
      npartsub+=f[ip][jp][kp];
    }
    double nden_s=npartsub/subboxsize;
    pshot[isub][jsub][ksub]=(double)1/nden_s;
    double norm_s=nden_s*sqrt(subboxsize);
    double pixden_s=npartsub/npsub[0]/npsub[1]/npsub[2];
    cout << pixden_s << endl;
    for (int i=0;i<npsub[0];i++) for (int j=0;j<npsub[1];j++) for (int k=0;k<npsub[2];k++) {
      int ip=i+isub*npsub[0], jp=j+jsub*npsub[1], kp=k+ksub*npsub[2];
      f[ip][jp][kp]=(f[ip][jp][kp]-pixden_s)/norm_s;
    }
  }
}

void smooth_field(double ***f, int *npix, const double& s) {
  fourier::westward_ho_3d(f,npix[0],npix[1],npix[2]);
  for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) for (int k=0;k<=npix[2]/2;k++) {
    int i1=i; if (i>npix[0]/2) i1=i-npix[0]; 
    int j1=j; if (j>npix[1]/2) j1=j-npix[1]; 
    double r=sqrt((double)((long)i1*i1+(long)j1*j1+(long)k*k))*2*pi, filter=gauss(r*s);   // (long)
    f[i][j][2*k]*=filter;
    f[i][j][2*k+1]*=filter;
  }
  fourier::eastward_ho_3d(f,npix[0],npix[1],npix[2]);
}

void smooth_field_sub(double ***f, int *npix, const double& s, const int& nsub) {
  int npsub[3]={npix[0]/nsub,npix[1]/nsub,npix[2]/nsub};
  for (int isub=0;isub<nsub;isub++) for (int jsub=0;jsub<nsub;jsub++) for (int ksub=0;ksub<nsub;ksub++) {
    double ***fsub=initmat3(npsub[0],npsub[1],npsub[2]+2);
    for (int i=0;i<npsub[0];i++) for (int j=0;j<npsub[1];j++) for (int k=0;k<npsub[2];k++) {
       fsub[i][j][k]=f[i+isub*npsub[0]][j+jsub*npsub[1]][k+ksub*npsub[2]];
    }
    fourier::westward_ho_3d(fsub,npsub[0],npsub[1],npsub[2]);
    for (int i=0;i<npsub[0];i++) for (int j=0;j<npsub[1];j++) for (int k=0;k<=npsub[2]/2;k++) {
      int i1=i; if (i>npsub[0]/2) i1=i-npsub[0]; 
      int j1=j; if (j>npsub[1]/2) j1=j-npsub[1]; 
      double r=sqrt((double)((long)i1*i1+(long)j1*j1+(long)k*k))*2*pi, filter=gauss(r*s);   // (long)
      fsub[i][j][2*k]*=filter;
      fsub[i][j][2*k+1]*=filter;
    }
    fourier::eastward_ho_3d(fsub,npsub[0],npsub[1],npsub[2]);
    for (int i=0;i<npsub[0];i++) for (int j=0;j<npsub[1];j++) for (int k=0;k<npsub[2];k++) {
      f[i+isub*npsub[0]][j+jsub*npsub[1]][k+ksub*npsub[2]]=fsub[i][j][k];
    }
    delete [] **fsub; delete [] *fsub; delete [] fsub;
  }
}

void deriv_field(double ***f, double ****v, int *npix, double *len) {
  fourier::westward_ho_3d(f,npix[0],npix[1],npix[2]);
  for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) for (int k=0;k<=npix[2]/2;k++) {
    int i1=i; if (i>npix[0]/2) i1=i-npix[0]; 
    int j1=j; if (j>npix[1]/2) j1=j-npix[1]; 
    double kv[3]={2*pi*i1/len[0],2*pi*j1/len[1],2*pi*k/len[2]}, kabs=norm(kv,3);
    if (fabs(kabs)>1e-10) {
      double kabs2=kabs*kabs;	
      for (int d=0;d<3;d++) {
	v[d][i][j][2*k]=kv[d]/kabs2*f[i][j][2*k+1];
	v[d][i][j][2*k+1]=-kv[d]/kabs2*f[i][j][2*k];
      }
    }
  }
  fourier::eastward_ho_3d(f,npix[0],npix[1],npix[2]);
  for (int d=0;d<3;d++) fourier::eastward_ho_3d(v[d],npix[0],npix[1],npix[2]);
  
  double sig[3]={0,0,0};
  for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) for (int k=0;k<npix[2];k++) for (int d=0;d<3;d++) sig[d]+=pow(v[d][i][j][k],2);
  cout << "sigv: ";
  for (int d=0;d<3;d++) cout << sqrt(sig[d]/(long)npix[0]/npix[1]/npix[2]) << " ";   // (long)
  cout << endl;
}

void deriv_field_sub(double ***f, double ****v, int *npix, double *lensub, const int& nsub) {
  int npsub[3]={npix[0]/nsub,npix[1]/nsub,npix[2]/nsub};
  for (int isub=0;isub<nsub;isub++) for (int jsub=0;jsub<nsub;jsub++) for (int ksub=0;ksub<nsub;ksub++) {
    double ***fsub=initmat3(npsub[0],npsub[1],npsub[2]+2);
    double ****vsub=initmat4(3,npsub[0],npsub[1],npsub[2]+2);
    for (int i=0;i<npsub[0];i++) for (int j=0;j<npsub[1];j++) for (int k=0;k<npsub[2];k++) {
       fsub[i][j][k]=f[i+isub*npsub[0]][j+jsub*npsub[1]][k+ksub*npsub[2]];
    }
    fourier::westward_ho_3d(fsub,npsub[0],npsub[1],npsub[2]);
    for (int i=0;i<npsub[0];i++) for (int j=0;j<npsub[1];j++) for (int k=0;k<=npsub[2]/2;k++) {
      int i1=i; if (i>npsub[0]/2) i1=i-npsub[0]; 
      int j1=j; if (j>npsub[1]/2) j1=j-npsub[1]; 
      double kv[3]={2*pi*i1/lensub[0],2*pi*j1/lensub[1],2*pi*k/lensub[2]}, kabs=norm(kv,3);
      if (fabs(kabs)>1e-10) {
	double kabs2=kabs*kabs;	
	for (int d=0;d<3;d++) {
	  vsub[d][i][j][2*k]=kv[d]/kabs2*fsub[i][j][2*k+1];
	  vsub[d][i][j][2*k+1]=-kv[d]/kabs2*fsub[i][j][2*k];
	}
      }
    }
    fourier::eastward_ho_3d(fsub,npsub[0],npsub[1],npsub[2]);
    for (int d=0;d<3;d++) fourier::eastward_ho_3d(vsub[d],npsub[0],npsub[1],npsub[2]);
    for (int i=0;i<npsub[0];i++) for (int j=0;j<npsub[1];j++) for (int k=0;k<npsub[2];k++) {
      f[i+isub*npsub[0]][j+jsub*npsub[1]][k+ksub*npsub[2]]=fsub[i][j][k];
      for (int d=0;d<3;d++) v[d][i+isub*npsub[0]][j+jsub*npsub[1]][k+ksub*npsub[2]]=vsub[d][i][j][k];
    }
    delete [] **fsub; delete [] *fsub; delete [] fsub;
    delete [] ***vsub; delete [] **vsub; delete [] *vsub; delete [] vsub;
  }
  
  double sig[3]={0,0,0};
  for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) for (int k=0;k<npix[2];k++) for (int d=0;d<3;d++) sig[d]+=pow(v[d][i][j][k],2);
  cout << "sigv: ";
  for (int d=0;d<3;d++) cout << sqrt(sig[d]/(long)npix[0]/npix[1]/npix[2]) << " ";   // (long)
  cout << endl;
}

void iso_shiftfield(double ****v, int *npix, double *len, const double& fz) {
  cout << "make shift field isotropic";
  for (int d=0;d<3;d++) fourier::westward_ho_3d(v[d],npix[0],npix[1],npix[2]);
  for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) for (int k=0;k<=npix[2]/2;k++) {
    int i1=i; if (i>npix[0]/2) i1=i-npix[0]; 
    int j1=j; if (j>npix[1]/2) j1=j-npix[1]; 
    double kv[3]={2*pi*i1/len[0],2*pi*j1/len[1],2*pi*k/len[2]}, kabs=norm(kv,3);
    if (fabs(kabs)>1e-10) {
      double mu=kv[2]/kabs;
      for (int d=0;d<3;d++) {
	v[d][i][j][2*k]/=(1+fz*mu*mu);
	v[d][i][j][2*k+1]/=(1+fz*mu*mu);
      }
    }
  }
  for (int d=0;d<3;d++) fourier::eastward_ho_3d(v[d],npix[0],npix[1],npix[2]);
}

void reconstruct_grid(double ***f, double ****v, int *npix, const double& eps, double *len) {
  int ip[3],ip1[3];
  double fp[3];
  
  for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) for (int k=0;k<npix[2];k++) f[i][j][k]=0.;
  for (int i=0;i<npix[0];i++) {
    for (int j=0;j<npix[1];j++) {
      for (int k=0;k<npix[2];k++) {
	fp[0]=i+eps*v[0][i][j][k]/len[0]*npix[0];
	fp[1]=j+eps*v[1][i][j][k]/len[1]*npix[1];
	fp[2]=k+eps*v[2][i][j][k]/len[2]*npix[2];
	for (int l=0;l<3;l++) {
	  if (fp[l]<0) fp[l]+=npix[l];
	  ip[l]=int(fp[l])%npix[l];
	  fp[l]-=int(fp[l]);
	  ip1[l]=(ip[l]+1)%npix[l];
	}
	f[ip[0]][ip[1]][ip[2]]+=(1-fp[0])*(1-fp[1])*(1-fp[2]);
	f[ip1[0]][ip[1]][ip[2]]+=fp[0]*(1-fp[1])*(1-fp[2]);
	f[ip[0]][ip1[1]][ip[2]]+=(1-fp[0])*fp[1]*(1-fp[2]);
	f[ip[0]][ip[1]][ip1[2]]+=(1-fp[0])*(1-fp[1])*fp[2];
	f[ip1[0]][ip1[1]][ip[2]]+=fp[0]*fp[1]*(1-fp[2]);
	f[ip[0]][ip1[1]][ip1[2]]+=(1-fp[0])*fp[1]*fp[2];
	f[ip1[0]][ip[1]][ip1[2]]+=fp[0]*(1-fp[1])*fp[2];
	f[ip1[0]][ip1[1]][ip1[2]]+=fp[0]*fp[1]*fp[2];
      }
    }
  }
}

void calc_pk(double ***f, double *len, int *npix, 
	     const double& kmin, const double& kmax, const int& nkbin, const double& shotn, char* outf) {
  
  fourier::westward_ho_3d(f,npix[0],npix[1],npix[2]);
  
  const int lmax=3;
  double nbin[nkbin], dnbin[nkbin], **dpk=initmat2(nkbin,lmax), **pk=initmat2(nkbin,lmax);
  for (int i=0;i<nkbin;i++) nbin[i]=0;
  double knyq=2*pi/len[0]*npix[0]/2;
  
  for (int k=0;k<=npix[2]/2;k++) {
    for (int i=0;i<nkbin;i++) {dnbin[i]=0; for (int l=0;l<lmax;l++) dpk[i][l]=0.;}
    for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) {
      int i1=i; if (i>npix[0]/2) i1=i-npix[0];
      int j1=j; if (j>npix[1]/2) j1=j-npix[1];
      double kv[3]={(float)i1/(npix[0]/2),(float)j1/(npix[1]/2),(float)k/(npix[2]/2)}, kabs=norm(kv,3);
      if (kabs>=kmin&&kabs<kmax) {
	double mu=kv[2]/kabs;
	int ibin=int((kabs-kmin)/(kmax-kmin)*nkbin);
	//int ibin=int(log(kabs/kmin)/log(kmax/kmin)*nkbin);
	double wkx,wky,wkz,wk2;
	if (i1*j1*k==0) wk2=1.; else {
	  wkx=sin(pi*i1/npix[0])/(pi*i1/npix[0]);
	  wky=sin(pi*j1/npix[1])/(pi*j1/npix[1]);
	  wkz=sin(pi*k/npix[2])/(pi*k/npix[2]);
	  wk2=pow(wkx*wky*wkz,4.0);
	  //wk2=1.;
	}
	double pl[3]={1,5./2.*(3*mu*mu-1),9./8.*(35*pow(mu,4)-30*mu*mu+3)};
	for (int l=0;l<lmax;l++) dpk[ibin][l]+=(pow(f[i][j][2*k],2)+pow(f[i][j][2*k+1],2))/wk2*pl[l];
	dnbin[ibin]++;
      }
    }
    int fac=1; if (k>0) fac=2;
    for (int i=0;i<nkbin;i++) {nbin[i]+=dnbin[i]*fac; for (int l=0;l<lmax;l++) pk[i][l]+=dpk[i][l]*fac;}
  }
  double x[nkbin];
  for (int i=0;i<nkbin;i++) {
    x[i]=knyq*(kmin+(kmax-kmin)/nkbin*(i+0.5));
    if (nbin[i]>0) {
      for (int l=0;l<lmax;l++) pk[i][l]/=nbin[i];
      pk[i][0]-=shotn;
    }
  }
  ofwritefunc_2cols(outf,x,pk,nkbin,lmax);

  //!!! Skip FFT to go back to real space
  //fourier::eastward_ho_3d(f,npix[0],npix[1],npix[2]);
}

void calc_pk_sub(double ***f, double *lensub, int *npix, 
		 const double& kmin, const double& kmax, const int& nkbin, double ***shotn, char* outf, const int& nsub) {
  
  int npsub[3]={npix[0]/nsub,npix[1]/nsub,npix[2]/nsub};
  for (int isub=0;isub<nsub;isub++) for (int jsub=0;jsub<nsub;jsub++) for (int ksub=0;ksub<nsub;ksub++) {
    //cout << shotn[isub][jsub][ksub] << endl;
    int subnum=isub*nsub*nsub+jsub*nsub+ksub;
    char outfsub[150];
    sprintf(outfsub,"%s%d",outf,subnum);
    cout << outfsub << endl;
    double ***fsub=initmat3(npsub[0],npsub[1],npsub[2]+2);
    for (int i=0;i<npsub[0];i++) for (int j=0;j<npsub[1];j++) for (int k=0;k<npsub[2];k++) {
       fsub[i][j][k]=f[i+isub*npsub[0]][j+jsub*npsub[1]][k+ksub*npsub[2]];
    }
    fourier::westward_ho_3d(fsub,npsub[0],npsub[1],npsub[2]);
  
    const int lmax=3;
    double nbin[nkbin], dnbin[nkbin], **dpk=initmat2(nkbin,lmax), **pk=initmat2(nkbin,lmax);
    for (int i=0;i<nkbin;i++) nbin[i]=0;
    double knyq=2*pi/lensub[0]*npsub[0]/2;
  
    for (int k=0;k<=npsub[2]/2;k++) {
      for (int i=0;i<nkbin;i++) {dnbin[i]=0; for (int l=0;l<lmax;l++) dpk[i][l]=0.;}
      for (int i=0;i<npsub[0];i++) for (int j=0;j<npsub[1];j++) {
	  int i1=i; if (i>npsub[0]/2) i1=i-npsub[0];
	  int j1=j; if (j>npsub[1]/2) j1=j-npsub[1];
	  double kv[3]={(float)i1/(npsub[0]/2),(float)j1/(npsub[1]/2),(float)k/(npsub[2]/2)}, kabs=norm(kv,3);
	  if (kabs>=kmin&&kabs<kmax) {
	    double mu=kv[2]/kabs;
	    int ibin=int((kabs-kmin)/(kmax-kmin)*nkbin);
	    //int ibin=int(log(kabs/kmin)/log(kmax/kmin)*nkbin);
	    double wkx,wky,wkz,wk2;
	    if (i1*j1*k==0) wk2=1.; else {
	      wkx=sin(pi*i1/npsub[0])/(pi*i1/npsub[0]);
	      wky=sin(pi*j1/npsub[1])/(pi*j1/npsub[1]);
	      wkz=sin(pi*k/npsub[2])/(pi*k/npsub[2]);
	      wk2=pow(wkx*wky*wkz,4.0);
	      //wk2=1.;
	    }
	    double pl[3]={1,5./2.*(3*mu*mu-1),9./8.*(35*pow(mu,4)-30*mu*mu+3)};
	    for (int l=0;l<lmax;l++) dpk[ibin][l]+=(pow(fsub[i][j][2*k],2)+pow(fsub[i][j][2*k+1],2))/wk2*pl[l];
	    dnbin[ibin]++;
	  }
	}
      int fac=1; if (k>0) fac=2;
      for (int i=0;i<nkbin;i++) {nbin[i]+=dnbin[i]*fac; for (int l=0;l<lmax;l++) pk[i][l]+=dpk[i][l]*fac;}
    }
    double x[nkbin];
    for (int i=0;i<nkbin;i++) {
      x[i]=knyq*(kmin+(kmax-kmin)/nkbin*(i+0.5));
      if (nbin[i]>0) {
	for (int l=0;l<lmax;l++) pk[i][l]/=nbin[i];
	pk[i][0]-=shotn[isub][jsub][ksub];
      }
    }
    ofwritefunc_2cols(outfsub,x,pk,nkbin,lmax);
    //!!! Skip FFT to go back to real space
    //fourier::eastward_ho_3d(fsub,npsub[0],npsub[1],npsub[2]);
    delete [] *dpk; delete [] dpk;
    delete [] *pk; delete [] pk;
    delete [] **fsub; delete [] *fsub; delete [] fsub;
  }
}

void calc_pk_cross(double ***f, double ***f2, double *len, int *npix, 
		   const double& kmin, const double& kmax, const int& nkbin, const double& shotn, const double& shotn2, char* outf) {
  
  fourier::westward_ho_3d(f,npix[0],npix[1],npix[2]);
  fourier::westward_ho_3d(f2,npix[0],npix[1],npix[2]);
  
  const int lmax=3;
  const int npow=3;
  double nbin[nkbin], dnbin[nkbin], ***dpk=initmat3(nkbin,lmax,npow), ***pk=initmat3(nkbin,lmax,npow);
  for (int i=0;i<nkbin;i++) nbin[i]=0;
  double knyq=2*pi/len[0]*npix[0]/2;
  
  for (int k=0;k<=npix[2]/2;k++) {
    for (int i=0;i<nkbin;i++) {dnbin[i]=0; for (int l=0;l<lmax;l++) for (int m=0;m<npow;m++) dpk[i][l][m]=0.;}
    for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) {
      int i1=i; if (i>npix[0]/2) i1=i-npix[0];
      int j1=j; if (j>npix[1]/2) j1=j-npix[1];
      double kv[3]={(float)i1/(npix[0]/2),(float)j1/(npix[1]/2),(float)k/(npix[2]/2)}, kabs=norm(kv,3);
      if (kabs>=kmin&&kabs<kmax) {
	double mu=kv[2]/kabs;
	int ibin=int((kabs-kmin)/(kmax-kmin)*nkbin);
	//int ibin=int(log(kabs/kmin)/log(kmax/kmin)*nkbin);
	double wkx,wky,wkz,wk2;
	if (i1*j1*k==0) wk2=1.; else {
	  wkx=sin(pi*i1/npix[0])/(pi*i1/npix[0]);
	  wky=sin(pi*j1/npix[1])/(pi*j1/npix[1]);
	  wkz=sin(pi*k/npix[2])/(pi*k/npix[2]);
	  wk2=pow(wkx*wky*wkz,4.0);
	  //wk2=1.;
	}
	double pl[3]={1,5./2.*(3*mu*mu-1),9./8.*(35*pow(mu,4)-30*mu*mu+3)};
	for (int l=0;l<lmax;l++) {
	  dpk[ibin][l][0]+=(pow(f[i][j][2*k],2)+pow(f[i][j][2*k+1],2))/wk2*pl[l];
	  dpk[ibin][l][1]+=(pow(f2[i][j][2*k],2)+pow(f2[i][j][2*k+1],2))/wk2*pl[l];
	  dpk[ibin][l][2]+=(f[i][j][2*k]*f2[i][j][2*k]+f[i][j][2*k+1]*f2[i][j][2*k+1])/wk2*pl[l];
	}
	dnbin[ibin]++;
      }
    }
    int fac=1; if (k>0) fac=2;
    for (int i=0;i<nkbin;i++) {nbin[i]+=dnbin[i]*fac; for (int l=0;l<lmax;l++) for (int m=0;m<npow;m++) pk[i][l][m]+=dpk[i][l][m]*fac;}
  }
  double x[nkbin];
  for (int i=0;i<nkbin;i++) {
    x[i]=knyq*(kmin+(kmax-kmin)/nkbin*(i+0.5));
    if (nbin[i]>0) {
      for (int l=0;l<lmax;l++) for (int m=0;m<npow;m++) pk[i][l][m]/=nbin[i];
      pk[i][0][0]-=shotn;
      pk[i][0][1]-=shotn2;
    }
  }
  ofwritefunc_3cols(outf,x,pk,nkbin,lmax,npow);

  //!!! Skip FFT to go back to real space
  //fourier::eastward_ho_3d(f,npix[0],npix[1],npix[2]);

  delete [] **dpk; delete [] *dpk; delete [] dpk;
  delete [] **pk; delete [] *pk; delete [] pk;

}

void calc_pk_cross_sub(double ***f, double ***f2, double *lensub, int *npix, 
		       const double& kmin, const double& kmax, const int& nkbin, double ***shotn, double ***shotn2, char* outf, const int& nsub) {
  
  int npsub[3]={npix[0]/nsub,npix[1]/nsub,npix[2]/nsub};
  for (int isub=0;isub<nsub;isub++) for (int jsub=0;jsub<nsub;jsub++) for (int ksub=0;ksub<nsub;ksub++) {
    //cout << shotn[isub][jsub][ksub] << endl;
    int subnum=isub*nsub*nsub+jsub*nsub+ksub;
    char outfsub[150];
    sprintf(outfsub,"%s%d",outf,subnum);
    cout << outfsub << endl;
    double ***fsub=initmat3(npsub[0],npsub[1],npsub[2]+2);
    double ***fsub2=initmat3(npsub[0],npsub[1],npsub[2]+2);
    for (int i=0;i<npsub[0];i++) for (int j=0;j<npsub[1];j++) for (int k=0;k<npsub[2];k++) {
       fsub[i][j][k]=f[i+isub*npsub[0]][j+jsub*npsub[1]][k+ksub*npsub[2]];
       fsub2[i][j][k]=f2[i+isub*npsub[0]][j+jsub*npsub[1]][k+ksub*npsub[2]];
    }
    fourier::westward_ho_3d(fsub,npsub[0],npsub[1],npsub[2]);
    fourier::westward_ho_3d(fsub2,npsub[0],npsub[1],npsub[2]);
  
    const int lmax=3;
    const int npow=3;
    double nbin[nkbin], dnbin[nkbin], ***dpk=initmat3(nkbin,lmax,npow), ***pk=initmat3(nkbin,lmax,npow);
    for (int i=0;i<nkbin;i++) nbin[i]=0;
    double knyq=2*pi/lensub[0]*npsub[0]/2;
  
    for (int k=0;k<=npsub[2]/2;k++) {
      for (int i=0;i<nkbin;i++) {dnbin[i]=0; for (int l=0;l<lmax;l++) for (int m=0;m<npow;m++) dpk[i][l][m]=0.;}
      for (int i=0;i<npsub[0];i++) for (int j=0;j<npsub[1];j++) {
	  int i1=i; if (i>npsub[0]/2) i1=i-npsub[0];
	  int j1=j; if (j>npsub[1]/2) j1=j-npsub[1];
	  double kv[3]={(float)i1/(npsub[0]/2),(float)j1/(npsub[1]/2),(float)k/(npsub[2]/2)}, kabs=norm(kv,3);
	  if (kabs>=kmin&&kabs<kmax) {
	    double mu=kv[2]/kabs;
	    int ibin=int((kabs-kmin)/(kmax-kmin)*nkbin);
	    //int ibin=int(log(kabs/kmin)/log(kmax/kmin)*nkbin);
	    double wkx,wky,wkz,wk2;
	    if (i1*j1*k==0) wk2=1.; else {
	      wkx=sin(pi*i1/npsub[0])/(pi*i1/npsub[0]);
	      wky=sin(pi*j1/npsub[1])/(pi*j1/npsub[1]);
	      wkz=sin(pi*k/npsub[2])/(pi*k/npsub[2]);
	      wk2=pow(wkx*wky*wkz,4.0);
	      //wk2=1.;
	    }
	    double pl[3]={1,5./2.*(3*mu*mu-1),9./8.*(35*pow(mu,4)-30*mu*mu+3)};
	    for (int l=0;l<lmax;l++) {
	      dpk[ibin][l][0]+=(pow(fsub[i][j][2*k],2)+pow(fsub[i][j][2*k+1],2))/wk2*pl[l];
	      dpk[ibin][l][1]+=(pow(fsub2[i][j][2*k],2)+pow(fsub2[i][j][2*k+1],2))/wk2*pl[l];
	      dpk[ibin][l][2]+=(fsub[i][j][2*k]*fsub2[i][j][2*k]+fsub[i][j][2*k+1]*fsub2[i][j][2*k+1])/wk2*pl[l];
	    }
	    dnbin[ibin]++;
	  }
	}
      int fac=1; if (k>0) fac=2;
      for (int i=0;i<nkbin;i++) {nbin[i]+=dnbin[i]*fac; for (int l=0;l<lmax;l++) for (int m=0;m<npow;m++) pk[i][l][m]+=dpk[i][l][m]*fac;}
    }
    double x[nkbin];
    for (int i=0;i<nkbin;i++) {
      x[i]=knyq*(kmin+(kmax-kmin)/nkbin*(i+0.5));
      if (nbin[i]>0) {
	for (int l=0;l<lmax;l++) for (int m=0;m<npow;m++) pk[i][l][m]/=nbin[i];
	pk[i][0][0]-=shotn[isub][jsub][ksub];
	pk[i][0][1]-=shotn2[isub][jsub][ksub];
      }
    }
    ofwritefunc_3cols(outfsub,x,pk,nkbin,lmax,npow);
    //!!! Skip FFT to go back to real space
    //fourier::eastward_ho_3d(fsub,npsub[0],npsub[1],npsub[2]);
    delete [] **dpk; delete [] *pk; delete [] dpk;
    delete [] **pk; delete [] *pk; delete [] pk;
    delete [] **fsub; delete [] *fsub; delete [] fsub;
  }
}

void ofwritefunc(char*string, double*x, double*y, const int& n) {
  ofstream ofs;
  ofs.open(string,ios::out); 
  ofs.precision(10);
  if (!ofs) throw ios::failure("Failed to open output file");
  for (int i=0;i<n;i++) {
    ofs << x[i] << " " << y[i] << endl;
  }
  ofs.close();
}

void ofwritefunc_2cols(char*string, double*x, double**y, const int& n, const int& m) {
  ofstream ofs;
  ofs.open(string,ios::out); 
  ofs.precision(10);
  if (!ofs) throw ios::failure("Failed to open output file");
  for (int i=0;i<n;i++) {
    ofs << x[i] << " ";
    for (int j=0;j<m;j++) ofs << y[i][j] << " ";
    ofs << endl;
  }
  ofs.close();
}

void ofwritefunc_3cols(char*string, double*x, double***y, const int& n, const int& m, const int& l) {
  ofstream ofs;
  ofs.open(string,ios::out); 
  ofs.precision(10);
  if (!ofs) throw ios::failure("Failed to open output file");
  for (int i=0;i<n;i++) {
    ofs << x[i] << " ";
    for (int k=0;k<l;k++) for (int j=0;j<m;j++) ofs << y[i][j][k] << " ";
    ofs << endl;
  }
  ofs.close();
}

double norm(double*x, const int& n) {
  double sm;
  sm=0;
  for (int i=0;i<n;i++) {
    sm+=x[i]*x[i];
  }
  return sqrt(sm);
}

double gauss(const double& k) {
  return exp(-pow(k,2)/2.);
}

double hrate(const Parlist par) {
  // assuming flat universe
  double ommz=par.omm*pow(1+par.z,3);
  double omdz=(1-par.omm)*pow(1+par.z,3*(1+par.w));
  return sqrt(ommz+omdz);
}

void assign_density(double **d, double *len, const int& n, double ***f, int *npix) {
  double dmin[3]; //datarange(d,n,len,dmin);
  for (int j=0;j<3;j++) dmin[j]=0;
  int ip[3],ip1[3];
  double fp[3];
  //!!!
  //for (int i=0;i<npix[0];i++) for (int j=0;j<npix[1];j++) for (int k=0;k<npix[2];k++) f[i][j][k]=0.;
  for (int i=0;i<n;i++) {
    for (int j=0;j<3;j++) {
      ip[j]=floor((d[i][j]-dmin[j])/len[j]*npix[j]-0.5);
      fp[j]=(d[i][j]-dmin[j])/len[j]*npix[j]-0.5-ip[j];
      ip[j]=(ip[j]+npix[j])%npix[j];
      ip1[j]=(ip[j]+1)%npix[j];
    }
    f[ip[0]][ip[1]][ip[2]]+=(1-fp[0])*(1-fp[1])*(1-fp[2]);
    f[ip1[0]][ip[1]][ip[2]]+=fp[0]*(1-fp[1])*(1-fp[2]);
    f[ip[0]][ip1[1]][ip[2]]+=(1-fp[0])*fp[1]*(1-fp[2]);
    f[ip[0]][ip[1]][ip1[2]]+=(1-fp[0])*(1-fp[1])*fp[2];
    f[ip1[0]][ip1[1]][ip[2]]+=fp[0]*fp[1]*(1-fp[2]);
    f[ip[0]][ip1[1]][ip1[2]]+=(1-fp[0])*fp[1]*fp[2];
    f[ip1[0]][ip[1]][ip1[2]]+=fp[0]*(1-fp[1])*fp[2];
    f[ip1[0]][ip1[1]][ip1[2]]+=fp[0]*fp[1]*fp[2];
  }
}

#undef pi
