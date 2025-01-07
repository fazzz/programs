
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "Simpson_integ.h"
#include "EF.h"

double k_B=1.98723e-3;
double T=300;

double oneD_Gaussian_base(double x,double nyu,double Sigma) {
  double N;
  double pi;

  pi=acos(-1.0);

  N=-0.5*(x-nyu)*(x-nyu)/Sigma/Sigma;
  N=exp(N);

  //  N=1.0/sqrt(2.0*pi)/Sigma*N;

  return N;
}

double Simpson_integ_oneD_Gaussian_base(int N, double minx,double maxx,
					double nyu,double Sigma) {
  int i;
  double delta;
  double dx,x,*y,S1,S2,S=0.0;

  if ((N%2)==1) {
    printf("error in Simpson Integrator\n"); return 0;
  }

  dx=(maxx-minx)/(double)N;

  if (dx<0.0) {
    printf("error in Simpson Integrator\n"); return 0;
  }

  y=(double *)emalloc(sizeof(double)*N);
  S1=0.0;
  S2=0.0;

  for (i=0;i<N;++i) {
    x=minx+dx*(double)i;
    y[i]=oneD_Gaussian_base(x,nyu,Sigma);
  }
  for (i=1;i<N-1;i+=2) {
    S1+=y[i];
  }
  for (i=2;i<N-2;i+=2) {
    S2+=y[i];
  }
  S=(y[0]+4.0*S1+2.0*S2+y[N-1])*dx/3.0;

  free(y);

  return S;
}

double Simpson_integ_C_oneD_Gaussian_base(int N, double minx,double maxx,
					  double k, double x0,double nyu,double Sigma,
					  int periodicflag, double periodicity) {
  int i;
  double delta;
  double dx,x,*y,S1,S2,S=0.0;

  if ((N%2)==1) {
    printf("error in Simpson Integrator\n"); return 0;
  }

  dx=(maxx-minx)/(double)N;

  if (dx<0.0) {
    printf("error in Simpson Integrator\n"); return 0;
  }

  y=(double *)emalloc(sizeof(double)*N);
  S1=0.0;
  S2=0.0;

  for (i=0;i<N;++i) {
    x=minx+dx*(double)i;
    delta=fabs(x-x0);
    if (periodicflag==ON) if (fabs(delta)>0.5*periodicity) delta=periodicity-delta;
    y[i]=exp(-0.5*k*delta*delta/k_B/T)*oneD_Gaussian_base(x,nyu,Sigma);
  }
  for (i=1;i<N-1;i+=2) {
    S1+=y[i];
  }
  for (i=2;i<N-2;i+=2) {
    S2+=y[i];
  }
  S=(y[0]+4.0*S1+2.0*S2+y[N-1])*dx/3.0;

  free(y);

  return S;
}

double Simpson_integ_C_P_GB(int N, double minx,double maxx,
			    int K, double *omega_k, 
			    double k, double x0, 
			    double *nyu,double *Sigma,
			    int periodicflag, double periodicity) {
  int i;
  double S=0.0;

  for (i=0;i<K;++i) 
    S=omega_k[i]*Simpson_integ_C_oneD_Gaussian_base(N,minx,maxx,k,x0,nyu[i],Sigma[i],periodicflag,periodicity);

  return S;
}

double C_func(double k, double x0, double x, int periodicflag, double periodicity){
  double y,delta;

  delta=fabs(x-x0);
  if (periodicflag==ON) if (fabs(delta)>0.5*periodicity) delta=periodicity-delta;

  y=exp(-0.5*k*delta*delta/k_B/T);

  return y;
}
