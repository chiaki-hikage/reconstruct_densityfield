//fourier.h
#ifndef fourier_h__
#define fourier_h__
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<math.h>
#include<fftw3.h>
using namespace std;

namespace fourier {
  void westward_ho_3d(double ***f, const int& ni, const int& nj, const int& nk);
  void eastward_ho_3d(double ***f, const int& ni, const int& nj, const int& nk);
  void eastward_ho_3d_float(float ***f, const int& ni, const int& nj, const int& nk);
};

namespace fourier {

  void westward_ho_3d(double ***f, const int& ni, const int& nj, const int& nk) {
    fftw_plan plan;
    //cout << "FFT start" << endl;
    plan=fftw_plan_dft_r2c_3d(ni,nj,nk,f[0][0],(fftw_complex*)f[0][0],FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    //cout << "FFT finish" << endl;
    return;
  }

  void eastward_ho_3d(double ***f, const int& ni, const int& nj, const int& nk) {
    int i,j,k;
    fftw_plan plan;
    plan=fftw_plan_dft_c2r_3d(ni,nj,nk,(fftw_complex*)f[0][0],f[0][0],FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    for(int i=0;i<ni;i++) for(int j=0;j<nj;j++) for(int k=0;k<nk+2;k++) f[i][j][k]/=ni*nj*nk;
    return;
  }

  void eastward_ho_3d_float(float ***f, const int& ni, const int& nj, const int& nk) {
    int i,j,k;
    fftw_plan plan;
    plan=fftw_plan_dft_c2r_3d(ni,nj,nk,(fftw_complex*)f[0][0],(double*)f[0][0],FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    for(int i=0;i<ni;i++) for(int j=0;j<nj;j++) for(int k=0;k<nk+2;k++) f[i][j][k]/=ni*nj*nk;
    return;
  }

}

#endif // fourier_h__
