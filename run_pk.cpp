#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<complex>
#include<cstdio>
#include<omp.h>
using namespace std;

//int nthreads=48;
int nthreads=8;

int main(int argc, char** argv)
{
  char buf[200],inf[200], outf[200];
  int i,nreal,ign,isnap;
  float zout;
  //const int lbox=4000,npix=2048,npart=4096, nfiles=2560, nsub=8;
  const int lbox=500,npix=512,npart=512, nfiles=48, nsub=1;

  nreal=atoi(argv[3]);  
  cout << nreal << endl;

  // use snapshot 1
  isnap=1;

  if (isnap==0) {
    zout=1.48;
  } else if (isnap==1) {
    zout=1.02;
  } else if (isnap==2) {
    zout=0.48;
  } else if (isnap==3) {
    zout=0;
  }

  cout << isnap << " " << zout << endl;

  #pragma omp parallel for private(buf,inf,outf) schedule(dynamic,1)
  for (i=0;i<nreal;i++) {
    //if(sign==0) sprintf(buf, "gunzip %s.%d.gz", argv[1], i);
    //if(sign==1) sprintf(buf, "gzip %s.%d", argv[1], i);

    // reconstruction (Rs=10)
    sprintf(inf,"%srun%05d/snapdir_%03d/snapshot_%03d",argv[1],i+1,isnap,isnap);
    sprintf(outf,"%spk_mm_snap%03d_rec_rs%d_run%05d",argv[2],isnap,10,i+1);
    sprintf(buf, "./pk_recon -l%d -p%d -i%s -o%s -fy -r%d -n%d -s%d -z%4.2f -tz -Tsym -S%d -c%s", lbox, npix, inf, outf, 10, npart, nfiles, zout, nsub,"y");
    cout << buf << endl;
    system(buf);

    // reconstruction (Rs=15)
    sprintf(inf,"%srun%05d/snapdir_%03d/snapshot_%03d",argv[1],i+1,isnap,isnap);
    sprintf(outf,"%spk_mm_snap%03d_rec_rs%d_run%05d",argv[2],isnap,15,i+1);
    sprintf(buf, "./pk_recon -l%d -p%d -i%s -o%s -fy -r%d -n%d -s%d -z%4.2f -tz -Tsym -S%d -c%s", lbox, npix, inf, outf, 15, npart, nfiles, zout, nsub,"y");
    cout << buf << endl;
    system(buf);

    // reconstruction (Rs=20)
    sprintf(inf,"%srun%05d/snapdir_%03d/snapshot_%03d",argv[1],i+1,isnap,isnap);
    sprintf(outf,"%spk_mm_snap%03d_rec_rs%d_run%05d",argv[2],isnap,20,i+1);
    sprintf(buf, "./pk_recon -l%d -p%d -i%s -o%s -fy -r%d -n%d -s%d -z%4.2f -tz -Tsym -S%d -c%s", lbox, npix, inf, outf, 20, npart, nfiles, zout, nsub,"y");
    cout << buf << endl;
    system(buf);
  }

  return 0;
}
