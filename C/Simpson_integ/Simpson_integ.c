
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "Simpson_integ.h"
#include "EF.h"

double Simpson_integ_C_oneD_Gaussian(int N, double minx,double maxx,
				     double k, double x0,double nyu,double Sigma,double pi) {
  int i;
  double dx,x,*y,S1,S2,S=0.0;

  if ((N%2)==1) {
    printf("error in Simpson Integrator\n"); return 0;
  }

  dx=(maxx-minx)/(double)N;

  if (dx<0.0) {
    printf("error in Simpson Integrator\n"); return 0;
  }

  y=(double *)gcemalloc(sizeof(double)*N);
  S1=0.0;
  S2=0.0;

  for (i=0;i<N;++i) {
    x=minx+dx*(double)i;
    y[i]=0.5*k*(x-x0)*(x-x0)*oneD_Gaussian(x,nyu,Sigma,pi);
  }
  for (i=1;i<N-1;i+=2) {
    S1+=y[i];
  }
  for (i=2;i<N-2;i+=2) {
    S2+=y[i];
  }
  S=(y[0]+4.0*S1+2.0*S2+y[N-1])*dx/3.0;

  return S;
}

double Simpson_integ_C_x_oneD_Gaussian(int N, double minx,double maxx,
				       double k, double x0,double nyu,double Sigma,double pi) {
  int i;
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
    y[i]=0.5*k*(x-x0)*(x-x0)*oneD_Gaussian(x,nyu,Sigma,pi)*x;
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

double Simpson_integ_C_x2_oneD_Gaussian(int N, double minx,double maxx,
					double k, double x0,double nyu,double Sigma,double pi) {
  int i;
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
    y[i]=0.5*k*(x-x0)*(x-x0)*oneD_Gaussian(x,nyu,Sigma,pi)*(x-nyu)*(x-nyu);
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

double Simpson_integ_oneD_Gaussian(int N, double minx,double maxx,double nyu,double Sigma,double pi) {
  int i;
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
    y[i]=oneD_Gaussian(x,nyu,Sigma,pi);
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

double oneD_Gaussian(double x,double nyu,double Sigma,double pi) {
  double N;

  N=-0.5*(x-nyu)*(x-nyu)/Sigma/Sigma;
  N=exp(N);

  N=1.0/sqrt(2.0*pi)/Sigma*N;

  return N;
}
