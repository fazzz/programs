
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "DCA.h"
#include "EF.h"
#include "LA.h"

double DCAb_backpass(double *fA2,double *fB1,
		     double *fA1,double *fB2,
		     double *PA21,double *PB12,
		     double *V,double Q,double *beta,double *gamma,double *W) {
  int i,j;
  double *WPA21,*WPA21fA1,*WPB12,*WPB12fB2;
  double *PA21fA1,*PB12fB2,*PA21fA1_PB12fB2_beta;

  WPA21=(double *)gcemalloc(sizeof(double)*6*6);
  WPA21fA1=(double *)gcemalloc(sizeof(double)*6);
  WPB12=(double *)gcemalloc(sizeof(double)*6*6);
  WPB12fB2=(double *)gcemalloc(sizeof(double)*6);

  LA_mmult(W,PA21,WPA21,6);
  mvmult(WPA21,fA1,PA21fA1,6);
  LA_mmult(W,PB12,WPB12,6);
  mvmult(WPB12,fB2,WPB12fB2,6);

  for (i=0;i<6;++i) {
    fB1[i]=WPA21fA1[i]-WPB12fB2[i]+gamma[i];
    fA2[i]=-fB1[i];
  }

  PA21fA1=(double *)gcemalloc(sizeof(double)*6);
  PB12fB2=(double *)gcemalloc(sizeof(double)*6);
  PA21fA1_PB12fB2_beta=(double *)gcemalloc(sizeof(double)*6);

  mvmult(PA21,fA1,PA21fA1,6);
  mvmult(PB12,fB2,PB12fB2,6);

  for (i=0;i<6;++i) PA21fA1_PB12fB2_beta[i]=PA21fA1[i]-PB12fB2[i]+beta[i];

  return 1.0/V[2*6+2]*(Q-V[2]*PA21fA1_PB12fB2_beta[2]);
}
