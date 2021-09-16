//recon.h
#ifndef RECON_H
#define RECON_H
#include<iostream>
#include<math.h>
#include<fstream>
#include<stdlib.h>
#include<string.h>
using namespace std;

class Parlist {

 public:
  
  Parlist(){}
  ~Parlist(){}

  double len;          // side lengh
  int npix;          // grid number at each side
  int nsnap;          // number of snapshots
  int npart;          // cubic root of particle number
  int nsub;          // if nsub>1, the box is divided into nsub^3 subboxes
  char *outf;         // output filename
  char *inf;          // filename of input image (ascii)
  char *rec_flag;      // reconstruct or not
  char *space;       // real space or redshift space
  char *rectype;       // when space='z', symmetric reconstructon or isotropic reconstruction 
  char *calc_crossP;   // calculate cross spectrum between pre-P and post-P or not [y or n]
  double fz;       // growth rate
  double kmin, kmax; // kmin, kmax in unit of knyq
  double z,ahz;       // redshift
  int nkbin;         // number of k binning
  double rscale, eps;          // reconstruction parameter
  double omm, h0, w;     // cosmological parameters
  double fz_boost;   // artificial increment of growth rate

  /* set default values for global variables */				   
  void defaultval(){
    len=500.; 
    npart=512;
    npix=512;
    nsnap=48;
    nsub=1;
    outf="output/pk_test.dat";
    inf="input/dm.dat";
    rec_flag="y";
    rscale=10.;
    eps=1.;
    omm=0.3156;
    h0=0.6727;
    w=-1;
    space="r";
    rectype="s";
    z=1;
    fz=1;
    fz_boost=1;
    calc_crossP="y";
  }
  
  /* random command line into the global parameter list */
  void argument(int argc,char*argv[]) {
    char c;
   
    while((char)EOF!=(c=getopt(argc,argv,"l:p:i:o:f:r:n:s:z:t:T:g:S:F:c:"))) 
      switch(c) {
      case 'l': len      = atof(optarg); break;
      case 'p': npix     = atoi(optarg); break;
      case 'i': inf      =      optarg ; break;
      case 'o': outf     =      optarg ; break;
      case 'f': rec_flag =       optarg; break;
      case 'r': rscale   = atof(optarg); break;
      case 'n': npart    = atoi(optarg); break;
      case 's': nsnap    = atoi(optarg); break;
      case 'z': z        = atof(optarg); break;
      case 't': space  = optarg; break;
      case 'T': rectype  = optarg; break;
      case 'g': fz  = atof(optarg); break;
      case 'S': nsub    = atoi(optarg); break;
      case 'F': fz_boost  = atof(optarg); break;
      case 'c': calc_crossP  = optarg; break;
      }
  }
  
};

#endif // RECON_H
