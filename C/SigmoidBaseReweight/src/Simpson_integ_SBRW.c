
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "Simpson_integ_SBRW_2.h"
#include "EF.h"

double Simpson_integ_oneD_SBRW_C_SB(int N, double minx,double maxx,double x_i,double x_k,double bk,double h) {
  int i;
  double delta,delta2;
  double dx,x,*y,S1,S2,S=0.0;
  double dpi;

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

    delta=fabs(x-x_i);  
    if (delta>0.5*dpi)  delta=dpi-delta;

    delta2=fabs(x-x_k); 
    if (delta2>0.5*dpi) delta2=dpi-delta2;

    //    y[i]=exp(-0.5*bk*delta*delta-0.5*delta2*delta2/h);
    y[i]=exp(-0.5*bk*delta*delta)*(1.0+exp(-delta/h));
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

double Simpson_integ_oneD_SBRW_SB(int N, double minx,double maxx,double x_k,double h) {
  int i;
  double delta;
  double dx,x,*y,S1,S2,S=0.0;
  double dpi;

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

    delta=fabs(x-x_k); 
    if (delta>0.5*dpi) delta=dpi-delta;

    //    y[i]=exp(-0.5*delta*delta/h);
    y[i]=1.0/exp(-delta/h);
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

