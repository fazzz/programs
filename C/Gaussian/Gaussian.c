


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "Gaussian.h"
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

double mixed_twod_Gaussian(double *x,double **nyu,double ***Sigma,double *pi_k,int K,double pi) {
  int i;
  double N=0.0;

  for (i=0;i<K;++i){
    N+=pi_k[i]*twoD_Gaussian(x,nyu[i],Sigma[i],pi);
  }

  return N;
}

void de_ln_twoD_Mixed_Gaussian(double *x,double *f,double **nyu,double ***Sigma,double *pi_k,int K,double pi) {
  int i,j,k;
  double det=0.0;
  double a,b,c,d;
  double denominator=0.0;

  f[0]=0.0;
  f[1]=0.0;

  for (i=0;i<K;++i){
    det=Sigma[i][0][0]*Sigma[i][1][1]-Sigma[i][0][1]*Sigma[i][1][0];
    a=1.0/det*Sigma[i][1][1];
    b=1.0/det*(-1.0*Sigma[i][0][1]);
    c=1.0/det*(-1.0*Sigma[i][1][0]);
    d=1.0/det*Sigma[i][0][0];

    denominator+=pi_k[i]*twoD_Gaussian(x,nyu[i],Sigma[i],pi);

    f[0]+=pi_k[i]*twoD_Gaussian(x,nyu[i],Sigma[i],pi)*0.5*(2.0*a*(x[0]-nyu[i][0])+(b+c)*(x[1]-nyu[i][1]));
    f[1]+=pi_k[i]*twoD_Gaussian(x,nyu[i],Sigma[i],pi)*0.5*(2.0*d*(x[1]-nyu[i][1])+(b+c)*(x[0]-nyu[i][1]));
  }

  f[0]=f[0]/denominator;
  f[1]=f[1]/denominator;

}

double *Create_twoD_GaussianMap(double minx,double maxx,double dx,
				double miny,double maxy,double dy,
				double *nyu,double **Sigma,double pi) {
  int i,j,k;
  int numpoint=0;
  double *x;
  double *prob;
  
  x=(double *)gcemalloc(sizeof(double)*2);
  prob=(double *)gcemalloc(sizeof(double)*1);

  x[0]=minx;
  for (i=0;x[0]<maxx;++i){
    x[0]=minx+dx*i;
    x[1]=miny;
    for (j=0;x[1]<maxy;++j){
      x[1]=miny+dy*j;
      prob[numpoint]=twoD_Gaussian(x,nyu,Sigma,pi);
      ++numpoint;
      prob=(double *)gcerealloc(prob,sizeof(double)*numpoint);
    }
  }

  return prob;
}

double *Create_mixed_twoD_GaussianMap(double minx,double maxx,double dx,
				      double miny,double maxy,double dy,
				      double **nyu,double ***Sigma,
				      double *pi_k,int K,double pi) {
  int i,j,k;
  int numpoint=0;
  double *x;
  double *prob;
  
  x=(double *)gcemalloc(sizeof(double)*2);
  prob=(double *)gcemalloc(sizeof(double)*1);

  x[0]=minx;
  for (i=0;x[0]<maxx;++i){
    x[0]=minx+dx*i;
    x[1]=miny;
    for (j=0;x[1]<maxy;++j){
      x[1]=miny+dy*j;
      prob[numpoint]=mixed_twod_Gaussian(x,nyu,Sigma,pi_k,K,pi);
      ++numpoint;
      prob=(double *)gcerealloc(prob,sizeof(double)*numpoint);
    }
  }

  return prob;
}
