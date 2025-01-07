#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "const.h"
#include "dPCA.h"

void calc_dPCA(/*double *dihed_traj,*/int time, int numdihed,
	      double covR[MAXNUMDIHED*MAXNUMDIHED], double w[MAXNUMDIHED])
{
  int i,j,k;

  char jobz='V';char uplo='U';
  double work[MAXNUMDIHED];
  long int lda=numdihed,lwork=3*numdihed,info;

  double var;

  calc_dcov(/*dihed_traj,*/time,numdihed);

  for (i=0;i<numdihed*numdihed;++i)
    covR[i]=cov[i];

  dsyev_(&jobz, &uplo, &lda, covR, &lda, w, work, &lwork, &info);
  if (info!=0){
    printf("error; eigen vector calculation!!");
    exit(1);
  }

}

void calc_dcov(/*double *dihed_traj,*/ int time, int numdihed){
  int i,j,jj,k,kk,l;

  cov=malloc(sizeof(double)*numdihed*numdihed);
  dihed_ave=malloc(sizeof(double)*numdihed);

  for (i=0;i<time;++i){
    for (j=0;j<numdihed;++j){
      dihed_ave[j]=(i*dihed_ave[j]+dihed_traj[i*numdihed+j])/(i+1);
    }
  }

  for (i=0;i<time;++i){
    for (j=0;j<numdihed;++j){
      for (jj=0;jj<numdihed;++jj){
	cov[j*numdihed+jj]=(i*cov[j*numdihed+jj]
			   +(dihed_traj[i*numdihed +j]-dihed_ave[j ])
			   *(dihed_traj[i*numdihed+jj]-dihed_ave[jj]))/(i+1);
      }
    }
  }

}
