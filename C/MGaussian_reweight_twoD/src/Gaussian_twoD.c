
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "Gaussian_twoD.h"
#include "EF.h"

double twoD_Gaussian(double *x,double *nyu,double **Sigma,double pi) {
  int i,j,k;

  double N,**InvSigma,detSigma=0.0,*vec;
  double *InvSigmaxvec;

  InvSigma=(double **)gcemalloc(sizeof(double *)*2);
  for (i=0;i<2;++i) InvSigma[i]=(double *)gcemalloc(sizeof(double)*2);
  InvSigmaxvec=(double *)gcemalloc(sizeof(double)*2);

  vec=(double *)gcemalloc(sizeof(double)*2);

  detSigma=Sigma[0][0]*Sigma[1][1]-Sigma[0][1]*Sigma[1][0];

  InvSigma[0][0]=Sigma[1][1]/detSigma;
  InvSigma[0][1]=-Sigma[1][0]/detSigma;
  InvSigma[1][0]=-Sigma[0][1]/detSigma;
  InvSigma[1][1]=Sigma[0][0]/detSigma;

  vec[0]=x[0]-nyu[0];
  vec[1]=x[1]-nyu[1];

  for (i=0;i<2;++i) {
    for (j=0;j<2;++j) {
      InvSigmaxvec[i]=0.0;
    }
  }

  
  for (i=0;i<2;++i) {
    for (j=0;j<2;++j) {
      InvSigmaxvec[i]+=InvSigma[i][j]*vec[j];
    }
  }

  N=-0.5*(vec[0]*InvSigmaxvec[0]+vec[1]*InvSigmaxvec[1]);
  N=exp(N);

  N=1.0/(2.0*pi)/sqrt(detSigma)*N;

  return N;
}
