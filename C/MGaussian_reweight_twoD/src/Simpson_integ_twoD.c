
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "Simpson_integ_twoD.h"
#include "Gaussian_twoD.h"
#include "EF.h"

double k_B=1.98723e-3;
double T=300;

double Simpson_integ_twoD_Gaussian(int N, double *minx,double *maxx,double *nyu,double **Sigma,double pi) {
  int i,j;
  double dx1,dx2,x[2],**y,v,S=0.0,*temp;

  if ((N%2)==1) { printf("error in Simpson Integrator\n"); return 0; }

  dx1=(maxx[0]-minx[0])/(double)N;  dx2=(maxx[1]-minx[1])/(double)N;

  if (dx1<0.0 || dx2<0.0) { printf("error in Simpson Integrator\n"); return 0; }

  y=(double **)emalloc(sizeof(double)*N);  for (i=0;i<N;++i) y[i]=(double *)emalloc(sizeof(double)*N);
  temp=(double *)emalloc(sizeof(double)*N);

  for (i=0;i<N;++i) {
    x[0]=minx[0]+dx1*(double)i;
    for (j=0;j<N;++j) {
      x[1]=minx[1]+dx2*(double)j;
      y[i][j]=twoD_Gaussian(x,nyu,Sigma,pi);
    }
  }

  for (i=1;i<N-1;i+=2) {
    v=y[i][0]+y[i][N-1];
    for (j=1;j<N-1;j+=2) {
      v+=(2.0*y[i][j]+4.0*y[i][j+1]);
    }
    temp[i]=v;
  }
  v=temp[0]+temp[N-1];

  for (i=2;i<N-2;i+=2) {
    v+=2.0*temp[i]+4.0*temp[i+1];
  }

  S=v*dx1*dx2/9.0;

  free(y);
  free(temp);

  return S;
}

double Simpson_integ_C_twoD_Gaussian(int N, double *minx,double *maxx,
				     double *k, double *x0,double *nyu,double **Sigma,double pi,
				     int periodicflag, double periodicity) {
  int i,j;
  double delta[2];
  double dx1,dx2,x[2],**y,S=0.0,v,*temp;

  if ((N%2)==1) {  printf("error in Simpson Integrator\n"); return 0;  }

  dx1=(maxx[0]-minx[0])/(double)N;  dx2=(maxx[1]-minx[1])/(double)N;

  if (dx1<0.0 || dx2<0.0) { printf("error in Simpson Integrator\n"); return 0; }

  y=(double **)emalloc(sizeof(double)*N);  for (i=0;i<N;++i) y[i]=(double *)emalloc(sizeof(double)*N);
  temp=(double *)emalloc(sizeof(double)*N);

  for (i=0;i<N;++i) {
    x[0]=minx[0]+dx1*(double)i;
    delta[0]=fabs(x[0]-x0[0]);
    if (periodicflag==ON) { if (delta[0]>0.5*periodicity) delta[0]=periodicity-delta[0]; }
    for (j=0;j<N;++j) {
      x[1]=minx[1]+dx2*(double)j;
      delta[1]=fabs(x[1]-x0[1]);
      if (periodicflag==ON) { if (delta[1]>0.5*periodicity) delta[1]=periodicity-delta[1]; }
      y[i][j]=exp(-0.5*k[0]*delta[0]*delta[0]/k_B/T)*exp(-0.5*k[1]*delta[1]*delta[1]/k_B/T)
	*twoD_Gaussian(x,nyu,Sigma,pi);
    }
  }

  for (i=1;i<N-1;i+=2) {
    v=y[i][0]+y[i][N-1];
    for (j=1;j<N-1;j+=2) {
      v+=(2.0*y[i][j]+4.0*y[i][j+1]);
    }
    temp[i]=v;
  }
  v=temp[0]+temp[N-1];

  for (i=2;i<N-2;i+=2) {
    v+=2.0*temp[i]+4.0*temp[i+1];
  }

  S=v*dx1*dx2/9.0;

  free(y);
  free(temp);

  return S;
}

double *Simpson_integ_C_x_uk_twoD_Gaussian(int N, double *minx,double *maxx,
					   double *k, double *x0,double *nyu,double **Sigma,double pi,
					   int periodicflag, double periodicity) {
  int i,j;
  double delta[2];
  double dx1,dx2,x[2],**y1,**y2,v1,v2,*temp1,*temp2,*S/*[2]={0.0,0.0}}*/;

  double delta2;

  S=(double *)gcemalloc(sizeof(double)*2);
  for (i=0;i<2;++i) S[i]=0.0;

  if ((N%2)==1) { printf("error in Simpson Integrator\n"); return 0; }

  dx1=(maxx[0]-minx[0])/(double)N;  dx2=(maxx[1]-minx[1])/(double)N;

  if (dx1<0.0 || dx2<0.0) { printf("error in Simpson Integrator\n"); return 0; }

  y1=(double **)emalloc(sizeof(double)*N);  for (i=0;i<N;++i) y1[i]=(double *)emalloc(sizeof(double)*N);
  temp1=(double *)emalloc(sizeof(double)*N);

  y2=(double **)emalloc(sizeof(double)*N);  for (i=0;i<N;++i) y2[i]=(double *)emalloc(sizeof(double)*N);
  temp2=(double *)emalloc(sizeof(double)*N);

  for (i=0;i<N;++i) {
    x[0]=minx[0]+dx1*(double)i;
    delta[0]=fabs(x[0]-x0[0]);
    if (periodicflag==ON) { if (delta[0]>0.5*periodicity) delta[0]=periodicity-delta[0]; }
    for (j=0;j<N;++j) {
      x[1]=minx[1]+dx2*(double)j;
      delta[1]=fabs(x[1]-x0[1]);
      if (periodicflag==ON) { if (delta[1]>0.5*periodicity) delta[1]=periodicity-delta[1]; }
      y1[i][j]=exp(-0.5*k[0]*delta[0]*delta[0]/k_B/T)*exp(-0.5*k[1]*delta[1]*delta[1]/k_B/T)
	*twoD_Gaussian(x,nyu,Sigma,pi)*(x[0]-nyu[0]);
      y2[i][j]=exp(-0.5*k[0]*delta[0]*delta[0]/k_B/T)*exp(-0.5*k[1]*delta[1]*delta[1]/k_B/T)
	*twoD_Gaussian(x,nyu,Sigma,pi)*(x[1]-nyu[1]);
    }
  }

  for (i=1;i<N-1;i+=2) {
    v1=y1[i][0]+y1[i][N-1];
    v2=y2[i][0]+y2[i][N-1];
    for (j=1;j<N-1;j+=2) {
      v1+=(2.0*y1[i][j]+4.0*y1[i][j+1]);
      v2+=(2.0*y2[i][j]+4.0*y2[i][j+1]);
    }
    temp1[i]=v1;
    temp2[i]=v2;
  }
  v1=temp1[0]+temp1[N-1];
  v2=temp2[0]+temp2[N-1];

  for (i=2;i<N-2;i+=2) {
    v1+=2.0*temp1[i]+4.0*temp1[i+1];
    v2+=2.0*temp2[i]+4.0*temp2[i+1];
  }

  S[0]=v1*dx1*dx2/9.0;
  S[1]=v2*dx1*dx2/9.0;

  free(y1);
  free(y2);

  return S;
}

double **Simpson_integ_C_x_uk2_twoD_Gaussian(int N, double *minx,double *maxx,
					     double *k, double *x0,double *nyu,double **Sigma,double pi,
					     int periodicflag, double periodicity) {
  int i,j;
  double delta[2],deltax_nyu[2];
  double dx1,dx2,x[2],**y1,**y2,**y3,**y4,*temp1,*temp2,*temp3,*temp4,**S;
  double v1,v2,v3,v4;

  S=(double **)gcemalloc(sizeof(double *)*2);
  for (i=0;i<2;++i) S[i]=(double *)gcemalloc(sizeof(double)*2);
  for (i=0;i<2;++i) for (j=0;j<2;++j) S[i][j]=0.0;

  if ((N%2)==1) { printf("error in Simpson Integrator\n"); return 0; }

  dx1=(maxx[0]-minx[0])/(double)N;  dx2=(maxx[1]-minx[1])/(double)N;

  if (dx1<0.0 || dx2<0.0) { printf("error in Simpson Integrator\n"); return 0; }

  y1=(double **)emalloc(sizeof(double)*N);  for (i=0;i<N;++i) y1[i]=(double *)emalloc(sizeof(double)*N);
  temp1=(double *)emalloc(sizeof(double)*N);

  y2=(double **)emalloc(sizeof(double)*N);  for (i=0;i<N;++i) y2[i]=(double *)emalloc(sizeof(double)*N);
  temp2=(double *)emalloc(sizeof(double)*N);

  y3=(double **)emalloc(sizeof(double)*N);  for (i=0;i<N;++i) y3[i]=(double *)emalloc(sizeof(double)*N);
  temp3=(double *)emalloc(sizeof(double)*N);

  y4=(double **)emalloc(sizeof(double)*N);  for (i=0;i<N;++i) y4[i]=(double *)emalloc(sizeof(double)*N);
  temp4=(double *)emalloc(sizeof(double)*N);

  for (i=0;i<N;++i) {
    x[0]=minx[0]+dx1*(double)i;
    delta[0]=fabs(x[0]-x0[0]);
    if (periodicflag==ON) { if (delta[0]>0.5*periodicity) delta[0]=periodicity-delta[0]; }
    for (j=0;j<N;++j) {
      x[1]=minx[1]+dx2*(double)j;
      delta[1]=fabs(x[1]-x0[1]);
      if (periodicflag==ON) { if (delta[1]>0.5*periodicity) delta[1]=periodicity-delta[1]; }

      y1[i][j]=exp(-0.5*k[0]*delta[0]*delta[0]/k_B/T)*exp(-0.5*k[1]*delta[1]*delta[1]/k_B/T)
	*twoD_Gaussian(x,nyu,Sigma,pi)*(deltax_nyu[0]*deltax_nyu[0]-Sigma[0][0]);

      y2[i][j]=exp(-0.5*k[0]*delta[0]*delta[0]/k_B/T)*exp(-0.5*k[1]*delta[1]*delta[1]/k_B/T)
	*twoD_Gaussian(x,nyu,Sigma,pi)*(deltax_nyu[0]*deltax_nyu[1]-Sigma[0][1]);

      y3[i][j]=exp(-0.5*k[0]*delta[0]*delta[0]/k_B/T)*exp(-0.5*k[1]*delta[1]*delta[1]/k_B/T)
	*twoD_Gaussian(x,nyu,Sigma,pi)*(deltax_nyu[1]*deltax_nyu[0]-Sigma[1][0]);

      y4[i][j]=exp(-0.5*k[0]*delta[0]*delta[0]/k_B/T)*exp(-0.5*k[1]*delta[1]*delta[1]/k_B/T)
	*twoD_Gaussian(x,nyu,Sigma,pi)*(deltax_nyu[1]*deltax_nyu[1]-Sigma[1][1]);
    }
  }

  for (i=1;i<N-1;i+=2) {
    v1=y1[i][0]+y1[i][N-1];
    v2=y2[i][0]+y2[i][N-1];
    v3=y3[i][0]+y3[i][N-1];
    v4=y4[i][0]+y4[i][N-1];
    for (j=1;j<N-1;j+=2) {
      v1+=(2.0*y1[i][j]+4.0*y1[i][j+1]);
      v2+=(2.0*y2[i][j]+4.0*y2[i][j+1]);
      v3+=(2.0*y3[i][j]+4.0*y3[i][j+1]);
      v4+=(2.0*y4[i][j]+4.0*y4[i][j+1]);
    }
    temp1[i]=v1;
    temp2[i]=v2;
    temp3[i]=v3;
    temp4[i]=v4;
  }
  v1=temp1[0]+temp1[N-1];
  v2=temp2[0]+temp2[N-1];
  v3=temp3[0]+temp3[N-1];
  v4=temp4[0]+temp4[N-1];

  for (i=2;i<N-2;i+=2) {
    v1+=2.0*temp1[i]+4.0*temp1[i+1];
    v2+=2.0*temp2[i]+4.0*temp2[i+1];
    v3+=2.0*temp3[i]+4.0*temp3[i+1];
    v4+=2.0*temp4[i]+4.0*temp4[i+1];
  }

  S[0][0]=v1*dx1*dx2/9.0;
  S[0][1]=v2*dx1*dx2/9.0;
  S[1][0]=v3*dx1*dx2/9.0;
  S[1][1]=v4*dx1*dx2/9.0;

  free(y1);
  free(y2);
  free(y3);
  free(y4);

  return S;
}
