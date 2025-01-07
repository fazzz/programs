
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABAb.h"
#include "EF.h"
#include "../LA/LA.h"

void ABAbs_trans_mainpass(int stref,int endref, CLTb *clt,
			 double* PA1_t,double* PA2_t,double* PA21_t,double* bA1_t ,double* bA2_t,
			 double* PA1,double* PA2,double* PA21,double* bA1 ,double* bA2){
  double* transMat;

  transMat=(double *)gcemalloc(sizeof(double)*6*6);
  ABAbs_mak_transMat(transMat,clt,stref,endref);

  ABAbs_trans_b(bA1_t,transMat,bA1);
  ABAbs_trans_b(bA2_t,transMat,bA2);

  ABAbs_trans_P(PA1_t,transMat,PA1);
  ABAbs_trans_P(PA2_t,transMat,PA2);
  ABAbs_trans_P(PA21_t,transMat,PA21);

}

void ABAbs_trans_backpass(){

}

void ABAbs_trans_b(double* b_trans,double *TMat,double *b){

  mvmult(TMat,b,b_trans,6);
}

void ABAbs_trans_P(double* P_trans,double* TMat,double* P){
  double *TMatt,*TMP;

  TMatt=(double *)gcemalloc(sizeof(double)*6*6);
  TMP=(double *)gcemalloc(sizeof(double)*6*6);

  mtrans(TMat,TMatt,6);
  LA_mmult(TMat,P,TMP,6);
  LA_mmult(TMP,TMatt,P_trans,6);
}

void ABAbs_mak_transMat(double* transMat,CLTb *clt,int nstart, int nend){
  int i,j,n;
  double *transMattemp,*temp;

  transMattemp=(double *)gcemalloc(sizeof(double)*6*6);
  temp=(double *)gcemalloc(sizeof(double)*6*6);

  for (i=0;i<6;++i)
    for (j=0;j<6;++j)
      transMat[i*6+j]=clt[nstart].TM[i][j];
      
  for (n=clt[nstart].nNumClutOfParent-1;n!=nend;n=clt[n].nNumClutOfParent-1) {
    for (i=0;i<6;++i)
      for (j=0;j<6;++j)
	temp[i*6+j]=clt[n].TM[i][j];
    mvmult(transMat,temp,transMattemp,6);
    for (i=0;i<6;++i)
      for (j=0;j<6;++j)
	transMat[i*6+j]=transMattemp[i*6+j];
  }
  /*********************************************/
  /* for (i=0;i<6;++i)			       */
  /*   for (j=0;j<6;++j)		       */
  /*     temp[i*6+j]=clt[n].TM[i][j]; */
  /* mvmult(transMat,temp,transMattemp,6);     */
  /* for (i=0;i<6;++i)			       */
  /*   for (j=0;j<6;++j)		       */
  /*     transMat[i*6+j]=transMattemp[i*6+j];  */
  /*********************************************/
}
