
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "DCA.h"
#include "EF.h"
#include "LA.h"

double DCAm_mainpass(double *PC1, double *PC2, 
		     double *PC12,double *PC21,
		     double *bC1, double *bC2,
		     double *PA1, double *PA2,
		     double *PA12,double *PA21,
		     double *bA1, double *bA2,
		     double *PB1, double *PB2,
		     double *PB12,double *PB21,
		     double *bB1, double *bB2,
		     double *Coacc, double Q) {
  double *V,*beta,*VS,*W,*gamma;

  V=(double *)gcemalloc(sizeof(double)*6*6);
  beta=(double *)gcemalloc(sizeof(double)*6);
  VS=(double *)gcemalloc(sizeof(double)*6);
  W=(double *)gcemalloc(sizeof(double)*6*6);
  gamma=(double *)gcemalloc(sizeof(double)*6);

  DCAm_cV(V,PA2,PB1);
  DCAm_cBeta(beta,bA2,bB1,Coacc);
  DCAm_cW(W,VS,V);
  DCAm_cGamma(gamma,beta,W,V,Q,VS);
  DCAm_cABI(PC1,PC2,PC12,PC21,bC1,bC2,
	    PA1,PA2,PA12,PA21,bA1,bA2,
	    PB1,PB2,PB12,PB21,PB1,bB2,
	    W,gamma);

}

double DCAm_cW(double *W,double *VS,double *V) {
  int i,j;
  double *STV;
  double *VSSTV;

  STV=(double *)gcemalloc(sizeof(double)*6);
  VSSTV=(double *)gcemalloc(sizeof(double)*6*6);

  for (i=0;i<6;++i) {
    VS[i]=V[i*6+2];
    STV[i]=V[2*6+i];
  }

  for (i=0;i<6;++i) for (j=0;j<6;++j) VSSTV[i*6+j]=VS[i]*STV[j]/V[2*6+2];
  for (i=0;i<6;++i) for (j=0;j<6;++j) W[i*6+j]=V[i*6+j]-VSSTV[i*6+j];

}

double DCAm_cGamma(double *gamma, double *beta, double *W,double *V, double Q, double *VS) {
  int i;
  double *Wbeta;

  Wbeta=(double *)gcemalloc(sizeof(double)*6);
  mvmult(W,beta,Wbeta,6);

  for (i=0;i<6;++i) gamma[i]=Wbeta[i]+VS[i]/V[2*6+2]*Q;

}

double DCAm_cV(double *V, double *PA, double *PB) {
  int i,j;
  double *Vtemp;

  Vtemp=(double *)gcemalloc(sizeof(double)*6*6);
  for (i=0;i<6;++i) for (j=0;j<6;++j) Vtemp[i*6+j]=PA[i*6+j]+PB[i*6+j];
  invm2(Vtemp,V,6);

}

double DCAm_cBeta(double *beta, double *b2A, double *b1B, double *Coacc) {
  int i;

  for (i=0;i<6;++i) beta[i]=b2A[i]-b1B[i]+Coacc[i];

}

double DCAm_cABI(double *PC1, double *PC2,
		 double *PC12,double *PC21,
		 double *bC1, double *bC2,
		 double *PA1, double *PA2,
		 double *PA12,double *PA21,
		 double *bA1, double *bA2,
		 double *PB1, double *PB2, 
		 double *PB12,double *PB21,
		 double *bB1, double *bB2,
		 double *W, double *gamma) {
  int i,j;
  double *PA12W,*PA12WP21,*PB21W,*PB21WP12;
  double *PA12gamma,*PB21gamma;

  PA12W=(double *)gcemalloc(sizeof(double)*6*6);
  PA12WP21=(double *)gcemalloc(sizeof(double)*6*6);
  PB21W=(double *)gcemalloc(sizeof(double)*6*6);
  PB21WP12=(double *)gcemalloc(sizeof(double)*6*6);

  PA12gamma=(double *)gcemalloc(sizeof(double)*6);
  PB21gamma=(double *)gcemalloc(sizeof(double)*6);
  
  LA_mmult(PA12,W,PA12W,6);
  //mtrans(PA12,PA21,6);
  LA_mmult(PA12W,PA21,PA12WP21,6);
  for (i=0;i<6;++i) for (j=0;j<6;++j) PC1[i*6+j]=PA1[i*6+j]-PA12WP21[i*6+j];
  
  //mtrans(PB12,PB21,6);
  LA_mmult(PB21W,PB12,PB21WP12,6);
  for (i=0;i<6;++i) for (j=0;j<6;++j) PC2[i*6+j]=PB1[i*6+j]-PB21WP12[i*6+j];

  LA_mmult(PB21W,PA21,PC21,6);
  mtrans(PC21,PC12,6);

  mvmult(PA12,gamma,PA12gamma,6);
  for (i=0;i<6;++i) bC1[i]=bA1[i]-PA12gamma[i];

  mvmult(PB21,gamma,PB21gamma,6);
  for (i=0;i<6;++i) bC2[i]=bB2[i]-PB21gamma[i];

}

