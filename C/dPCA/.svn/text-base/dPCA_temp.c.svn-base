#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "f2c.h"
#include "blaswrap.h"
#include "clapack.h"

#include "dPCA.h"

int calc_dPC(int numdihed, int numstep, int numpc){
  int i,j,jj,k,kk;
  double pi,sum;
  double *q,*cov,*ave;

  char jobz='V';char uplo='U';
  long int n,lda,lwork,info;
  double covR[2*MAXNUMDIHED*2*MAXNUMDIHED];
  double w[2*MAXNUMDIHED];
  double work[3*2*MAXNUMDIHED];
  FILE *outputfile;

  pi=acos(-1);
  n=2*numdihed;lda=n;lwork=3*n;

  if((q=malloc(sizeof(double)*numstep*numdihed*2))==NULL){
    printf("error:cannot allocate q");
    exit(1);
  }
  if((ave=malloc(sizeof(double)*numdihed*2))==NULL){
    printf("error:cannot allocate ave");
    exit(1);
  }
  if((cov=malloc(sizeof(double)*numdihed*2*numdihed*2))==NULL) {
    printf("error:cannot allocate cov");
    exit(1);
  }
  if((dPC=malloc(sizeof(double)*numstep*numpc))==NULL) {
    printf("error:cannot allocate cov");
    exit(1);
  }
  if((dPV=malloc(sizeof(double)*numdihed*2))==NULL) {
    printf("error:cannot allocate cov");
    exit(1);
  }

  ave[0]=0.0;ave[1]=0.0;
  for (i=0;i<numstep;++i) {
    for (j=0;j<numdihed;++j) {
      q[i*2*numdihed+j*2]=cos(dihedtraj[i*numdihed+j]/180.0*pi);
      q[i*2*numdihed+j*2+1]=sin(dihedtraj[i*numdihed+j]/180.0*pi);
      ave[j*2]=(i*ave[j*2]+q[i*2*numdihed+j*2])/(i+1);
      ave[j*2+1]=(i*ave[j*2+1]+q[i*2*numdihed+j*2+1])/(i+1);
    }
  }

  for (j=0;j<numdihed*2;++j)
    for (jj=0;jj<numdihed*2;++jj)
      cov[j*numdihed*2+jj]=0.0;
  for (i=0;i<numstep;++i){
    for (j=0;j<numdihed*2;++j){
      for (jj=0;jj<numdihed*2;++jj){
	cov[j*numdihed*2+jj]=(i*cov[j*numdihed*2+jj]+(q[i*numdihed*2+j]-ave[j])*(q[i*numdihed*2+jj]-ave[jj]))/(i+1);
      }
    }
  }
  
  for (i=0;i<2*numdihed*2*numdihed;++i)
    covR[i]=cov[i];

  dsyev_(&jobz,&uplo,&n,covR,&lda,w,work,&lwork,&info);
  if (info!=0){
    printf("error; eigen vector calculation!!");
    exit(1);
  }
  sum=0.0;
  for (i=0;i<numdihed*2;++i) {
    dPV[i]=w[numdihed*2-1-i];
    sum+=dPV[i];
  }
  for (i=0;i<numdihed*2;++i) {
    printf("%8.3lf %8.3lf\n",dPV[i],dPV[i]/sum*100.0);
  }  
  for (i=0;i<numstep;++i) {
    for (j=0;j<numpc;++j) {
      for (k=0;k<numdihed*2;++k) {
	dPC[i*numpc+j]=covR[numdihed*2*(numdihed*2-j-1)+k]*q[i*numdihed*2+k];
      }  
    }
  }

  free(ave);
  free(cov);

  return 1;
}

int scandtraj(FILE *inputfile,/*double *dihedtraj,*/int numdihed, int numstep) {
  int i,j;
  double d;

  if ((dihedtraj = malloc(sizeof(double)*numdihed*numstep))==NULL){
    printf("error:cannot allocate traj\n");
    exit(1);
  }

  for (i=0;i<numstep;++i) {
    for (j=0;j<numdihed;++j) {
      fscanf(inputfile,"%lf",&d);
      dihedtraj[i*numdihed+j] = d;
    }
  }

  return 1;
}
