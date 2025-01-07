
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "Simpson_integ_TPBRW.h"
#include "EF.h"

double Simpson_integ_oneD_TPBRW_C(int N, double minx,double maxx,double x_i,int K,double bk,double beta,
				  double a0, double *a_k, double *b_k) {
  int i,j;
  double delta,delta2;
  double dx,x,*y,S1,S2,S=0.0;
  double dpi;

  double k_x,f;

  dpi=2.0*acos(-1.0);

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

    f=0.5*a0;
    for (j=0;j<K;++j) {
      k_x=(double)(j+1)*x;
      while (k_x > dpi) k_x-=dpi;
      while (k_x <= 0.0) k_x+=dpi;
      f+=a_k[j]*sin(k_x)+b_k[j]*cos(k_x);
    }
    //    y[i]=exp(-1.0*beta*f);
    delta=fabs(x-x_i);  
    if (delta>0.5*dpi)  delta=dpi-delta;
    //    y[i]=y[i]*exp(-0.5*bk*delta*delta);
    y[i]=exp(-1.0*beta*f)*exp(-0.5*bk*delta*delta);
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

double Simpson_integ_oneD_TPBRW_C_sin(int N, double minx,double maxx,double x_i,int k,int K,double bk,double beta,
				      double a0, double *a_k, double *b_k) {
  int i,j;
  double delta,delta2;
  double dx,x,*y,S1,S2,S=0.0;
  double dpi;

  double k_x,f;

  dpi=2.0*acos(-1.0);

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

    f=0.5*a0;
    for (j=0;j<K;++j) {
      k_x=(double)(j+1)*x;
      while (k_x > dpi) k_x-=dpi;
      while (k_x <= 0.0) k_x+=dpi;
      f+=a_k[j]*sin(k_x)+b_k[j]*cos(k_x);
    }
    //    y[i]=exp(-1.0*beta*f);
    delta=fabs(x-x_i);  
    if (delta>0.5*dpi)  delta=dpi-delta;
    k_x=(double)(k+1)*x;
    while (k_x > dpi) k_x-=dpi;
    while (k_x <= 0.0) k_x+=dpi;
    //    y[i]=y[i]*sin(k_x)*exp(-0.5*bk*delta*delta);
    y[i]=sin(k_x)*exp(-1.0*beta*f)*exp(-0.5*bk*delta*delta);
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

double Simpson_integ_oneD_TPBRW_C_cos(int N, double minx,double maxx,double x_i,int k,int K,double bk,double beta,
				      double a0, double *a_k, double *b_k) {

  int i,j;
  double delta,delta2;
  double dx,x,*y,S1,S2,S=0.0;
  double dpi;

  double k_x,f;

  dpi=2.0*acos(-1.0);

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

    f=0.5*a0;
    for (j=0;j<K;++j) {
      k_x=(double)(j+1)*x;
      while (k_x > dpi) k_x-=dpi;
      while (k_x <= 0.0) k_x+=dpi;
      f+=a_k[j]*sin(k_x)+b_k[j]*cos(k_x);
    }
    //    y[i]=exp(-1.0*beta*f);
    delta=fabs(x-x_i);  
    if (delta>0.5*dpi)  delta=dpi-delta;
    k_x=(double)(k+1)*x;
    while (k_x > dpi) k_x-=dpi;
    while (k_x <= 0.0) k_x+=dpi;
    //    y[i]=y[i]*cos(k_x)*exp(-0.5*bk*delta*delta);
    y[i]=cos(k_x)*exp(-1.0*beta*f)*exp(-0.5*bk*delta*delta);
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

