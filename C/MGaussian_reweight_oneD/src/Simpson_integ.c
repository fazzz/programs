
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "Simpson_integ.h"
#include "EF.h"

double k_B=1.98723e-3;
double T=300;

double Simpson_integ_C_oneD_Gaussian(int N, double minx,double maxx,
				     double k, double x0,double nyu,double Sigma,double pi,
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
    if (periodicflag==ON) if (fabs(delta)>0.5*periodicity) delta=/*2.0**/periodicity-delta;
    //    y[i]=0.5*k*delta*delta*oneD_Gaussian(x,nyu,Sigma,pi);
    y[i]=exp(-0.5*k*delta*delta/k_B/T)*oneD_Gaussian(x,nyu,Sigma,pi);
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

double Simpson_integ_C_x_oneD_Gaussian(int N, double minx,double maxx,
				       double k, double x0,double nyu,double Sigma,double pi,
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
    if (periodicflag==ON) if (fabs(delta)>0.5*periodicity) delta=/*2.0**/periodicity-delta;
    //    y[i]=0.5*k*delta*delta*oneD_Gaussian(x,nyu,Sigma,pi)*x;
    y[i]=exp(-0.5*k*delta*delta/k_B/T)*oneD_Gaussian(x,nyu,Sigma,pi)*x;
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
					double k, double x0,double nyu,double Sigma,double pi,
					int periodicflag, double periodicity) {
  int i;
  double delta,deltax_nyu;
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
    if (periodicflag==ON) if (fabs(delta)>0.5*periodicity) delta=/*2.0**/periodicity-delta;
    deltax_nyu=fabs(x-nyu);
    /********************************************************************************************************/
    /* if (periodicflag==ON) if (fabs(deltax_nyu)>0.5*periodicity) deltax_nyu=periodicity-fabs(deltax_nyu); */
    /* //    y[i]=0.5*k*delta*delta*oneD_Gaussian(x,nyu,Sigma,pi)*deltax_nyu*deltax_nyu;		    */
    /* y[i]=exp(-0.5*k*delta*delta/k_B/T)*oneD_Gaussian(x,nyu,Sigma,pi)*deltax_nyu*deltax_nyu;		    */
    /********************************************************************************************************/
    y[i]=exp(-0.5*k*delta*delta/k_B/T)*oneD_Gaussian(x,nyu,Sigma,pi)*(x-nyu)*(x-nyu);
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

double Simpson_integ_C_x_uk_oneD_Gaussian(int N, double minx,double maxx,
					  double k, double x0,double nyu,double Sigma,double pi,
					  int periodicflag, double periodicity) {
  int i;
  double delta;
  double dx,x,*y,S1,S2,S=0.0;

  double delta2;

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
    if (periodicflag==ON) if (fabs(delta)>0.5*periodicity) delta=/*2.0**/periodicity-delta;
    //    y[i]=0.5*k*delta*delta*oneD_Gaussian(x,nyu,Sigma,pi)*x;
    /******************************************************************************************************/
    /* delta2=fabs(x-nyu);										  */
    /* if (periodicflag==ON) if (fabs(delta2)>0.5*periodicity) delta2=/\*2.0**\/periodicity-fabs(delta2); */
    /* y[i]=exp(-0.5*k*delta*delta/k_B/T)*oneD_Gaussian(x,nyu,Sigma,pi)*delta2;				  */
    /******************************************************************************************************/
    y[i]=exp(-0.5*k*delta*delta/k_B/T)*oneD_Gaussian(x,nyu,Sigma,pi)*(x-nyu);
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

double Simpson_integ_C_x_uk2_oneD_Gaussian(int N, double minx,double maxx,
					   double k, double x0,double nyu,double Sigma,double pi,
					   int periodicflag, double periodicity) {
  int i;
  double delta,deltax_nyu;
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
    if (periodicflag==ON) if (fabs(delta)>0.5*periodicity) delta=/*2.0**/periodicity-delta;
    //    deltax_nyu=x-nyu;
    deltax_nyu=fabs(x-nyu);
    //    if (periodicflag==ON) if (deltax_nyu>0.5*periodicity) deltax_nyu=/*2.0**/periodicity-deltax_nyu;
    //    y[i]=exp(-0.5*k*delta*delta/k_B/T)*oneD_Gaussian(x,nyu,Sigma,pi)*(deltax_nyu*deltax_nyu-Sigma*Sigma);
    //    y[i]=0.5*k*delta*delta*oneD_Gaussian(x,nyu,Sigma,pi)*deltax_nyu*deltax_nyu;
    //    if (periodicflag==ON) if (fabs(deltax_nyu)>0.5*periodicity) deltax_nyu=/*2.0**/periodicity-deltax_nyu;
    y[i]=exp(-0.5*k*delta*delta/k_B/T)*oneD_Gaussian(x,nyu,Sigma,pi)*(deltax_nyu*deltax_nyu-Sigma*Sigma);
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
