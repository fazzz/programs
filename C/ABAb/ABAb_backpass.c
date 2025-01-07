#include <stdio.h>
#include <math.h>

#include "ABAb.h"

void ABAbm_backpass(double *qacc,ABIb* abi,CLTb *clt,int numclut,double *qvel,double q_NVT,double qvel_NVT,double *qacc_NVT,double s_NVT,double KE,double KEobj,int MODE) {
  int i,j,k;
  int nparent;
  double UNIT=418.4070;

  for(i=0;i<6;++i) clt[0].Spacc[i]=0.0;
  for (i=1;i<numclut;++i) {
    nparent=clt[i].nNumClutOfParent-1;

    for(j=0;j<6;++j) {
      clt[i].Spacc[j] = 0.0;
      clt[i].preSpacc[j] = 0.0;
    }

    for (j=0;j<6;++j) 
      for (k=0;k<6;++k) 
	clt[i].preSpacc[j]+=clt[i].TM[k][j]*clt[nparent].Spacc[k];
    
    qacc[i] = 0.0;
    for(j=0;j<6;++j)
      qacc[i]-=abi[i].KG[j]*clt[i].preSpacc[j];
    qacc[i]+=abi[i].nyu;
    
    for(j=0;j<6;++j) clt[i].Spacc[j]+=clt[i].preSpacc[j]+clt[i].Coacc[j];
    clt[i].Spacc[2]+=qacc[i];
  }

  if (MODE == NVT) {
    for(i=0;i<numclut;++i) qacc[i] -= qvel_NVT/q_NVT*qvel[i];
    *qacc_NVT = 2.0/s_NVT*q_NVT*KE*UNIT-KEobj/s_NVT*q_NVT*UNIT+qvel_NVT*qvel_NVT/q_NVT;
  }
}

void ABAbm_backpass_TermOn(double *qacc,ABIb* abi,CLTb *clt,int numclut,double *qvel,double q_NVT,double qvel_NVT,double *qacc_NVT,double s_NVT,double *acc_Term,double *acc_Term2,double KE,double KEobj,int MODE) {
  int i,j,k;
  int nparent;
  double UNIT=418.4070;

  for(i=0;i<6;++i) clt[0].Spacc[i]=acc_Term2[i];
  for (i=1;i<numclut;++i) {
    nparent=clt[i].nNumClutOfParent-1;

    for(j=0;j<6;++j) {
      clt[i].Spacc[j] = 0.0;
      clt[i].preSpacc[j] = 0.0;
    }

    for (j=0;j<6;++j) 
      for (k=0;k<6;++k) 
	clt[i].preSpacc[j]+=clt[i].TM[k][j]*clt[nparent].Spacc[k];
    
    qacc[i] = 0.0;
    for(j=0;j<6;++j)
      qacc[i]-=abi[i].KG[j]*clt[i].preSpacc[j];
    qacc[i]+=abi[i].nyu;
    
    for(j=0;j<6;++j) clt[i].Spacc[j]+=clt[i].preSpacc[j]+clt[i].Coacc[j];
    clt[i].Spacc[2]+=qacc[i];
  }

  if (MODE == NVT) {
    for(i=0;i<numclut;++i) qacc[i] -= qvel_NVT/q_NVT*qvel[i];
    *qacc_NVT = 2.0/s_NVT*q_NVT*KE*UNIT-KEobj/s_NVT*q_NVT*UNIT+qvel_NVT*qvel_NVT/q_NVT;
  }
}

void ABAbm_backpass_TERM(double acc_Term2[6],double acc_Term[6],double vel_Term[6],double s_NVT,double svel_NVT,int MODE) {
  /******************************************************************/
  /* int i;							    */
  /* double add[6];						    */
  /* 								    */
  /* add[0]=0.0;						    */
  /* add[1]=0.0;						    */
  /* add[2]=0.0;						    */
  /* add[3]=(-vel_Term[2]*vel_Term[3+1]+vel_Term[1]*vel_Term[3+2]); */
  /* add[4]=( vel_Term[2]*vel_Term[3+0]-vel_Term[0]*vel_Term[3+2]); */
  /* add[5]=(-vel_Term[1]*vel_Term[3+0]+vel_Term[0]*vel_Term[3+1]); */
  /* 								    */
  /* if (MODE==NVT) {						    */
  /*   for (i=0;i<6;++i) {					    */
  /*     acc_Term2[i]-=dot_s_NVT/s_NVT*vel_Term[i];		    */
  /*     acc_Term[i]=acc_Term2[i]-add[i];			    */
  /*   }							    */
  /* }								    */
  /******************************************************************/
}
