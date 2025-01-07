
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABAb.h"
#include "EF.h"

double solverABA(double *qacc,double *qvel,CLT *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double q_NVT,double qvel_NVT,double *qacc_NVT,double s_NVT,double KE,double KEobj,int MODE){
  ABI *abi;

  abi=(ABI *)gcemalloc(sizeof(ABI)*numclt);

  ABAp_prepass(clt,qvel,numclt,numatom,crd);

  ABAs_forc(clt,frc,numclt);

  ABAm_mainpass(abi,clt,Q,numclt);
  ABAm_backpass(qacc,abi,clt,numclt,qvel,q_NVT,qvel_NVT,qacc_NVT,s_NVT,KE,KEobj,MODE);

}

double solverABA_TermOn(double *qacc,double *qvel,CLT *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double q_NVT,double qvel_NVT,double *qacc_NVT,double s_NVT,double *acc_Term,double *acc_Term2,double *vel_Term,double KE,double KEobj,int MODE){
  ABI *abi;

  abi=(ABI *)gcemalloc(sizeof(ABI)*numclt);

  ABAp_prepass_TermOn(clt,qvel,numclt,numatom,crd,vel_Term);

  ABAs_forc(clt,frc,numclt);

  ABAm_mainpass_TermOn(abi,clt,Q,numclt,acc_Term,acc_Term2,vel_Term);
  ABAm_backpass_TermOn(qacc,abi,clt,numclt,qvel,q_NVT,qvel_NVT,qacc_NVT,s_NVT,acc_Term,acc_Term2,KE,KEobj,MODE);

}

double solverABA_NH(double *qacc,double *qvel,CLT *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double s,double s_vel,double *s_acc,double tau2,double Temp,double TempB){
  ABI *abi;

  abi=(ABI *)gcemalloc(sizeof(ABI)*numclt);

  ABAp_prepass(clt,qvel,numclt,numatom,crd);

  ABAs_forc(clt,frc,numclt);

  ABAm_mainpass(abi,clt,Q,numclt);

  ABAm_backpass_NH(qacc,abi,clt,numclt,qvel,s,s_vel,s_acc,tau2,Temp,TempB);

}

double solverABA_TermOn_NH(double *qacc,double *qvel,CLT *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double s,double s_vel,double *s_acc,double tau2,double *acc_Term,double *acc_Term2,double *vel_Term,double Temp,double TempB){
  ABI *abi;

  abi=(ABI *)gcemalloc(sizeof(ABI)*numclt);

  ABAp_prepass_TermOn(clt,qvel,numclt,numatom,crd,vel_Term);

  ABAs_forc(clt,frc,numclt);

  ABAm_mainpass_TermOn_NH(abi,clt,Q,numclt,acc_Term,acc_Term2,vel_Term,s,s_vel);

  ABAm_backpass_TermOn_NH(qacc,abi,clt,numclt,qvel,s,s_vel,s_acc,tau2,acc_Term,acc_Term2,vel_Term,Temp,TempB);

}

double solverABA_NH_new(double *qacc,double *qvel,CLT *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double gzi,double *gzi_vel,double s,double *s_vel,double tau2,double Temp,double TempB){
  ABI *abi;

  abi=(ABI *)gcemalloc(sizeof(ABI)*numclt);

  ABAp_prepass(clt,qvel,numclt,numatom,crd);

  ABAs_forc(clt,frc,numclt);

  ABAm_mainpass(abi,clt,Q,numclt);

  ABAm_backpass_NH_new(qacc,abi,clt,numclt,qvel,gzi,gzi_vel,tau2,Temp,TempB);

}

double solverABA_TermOn_NH_new(double *qacc,double *qvel,CLT *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double gzi,double *gzi_vel,double s,double *s_vel,double tau2,double *acc_Term,double *acc_Term2,double *vel_Term,double Temp,double TempB){
  ABI *abi;

  abi=(ABI *)gcemalloc(sizeof(ABI)*numclt);

  ABAp_prepass_TermOn(clt,qvel,numclt,numatom,crd,vel_Term);

  ABAs_forc(clt,frc,numclt);

  ABAm_mainpass_TermOn_NH_new(abi,clt,Q,numclt,acc_Term,acc_Term2,vel_Term);

  ABAm_backpass_TermOn_NH_new(qacc,abi,clt,numclt,qvel,gzi,gzi_vel,tau2,acc_Term,acc_Term2,vel_Term,Temp,TempB);

}

double solverABA_NH_new_mvV(double *qacc,double *qvel,CLT *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double gzi){
  ABI *abi;

  abi=(ABI *)gcemalloc(sizeof(ABI)*numclt);

  ABAp_prepass(clt,qvel,numclt,numatom,crd);

  ABAs_forc(clt,frc,numclt);

  ABAm_mainpass(abi,clt,Q,numclt);

  ABAm_backpass_NH_new_mvV(qacc,abi,clt,numclt,qvel,gzi);

}

double solverABA_TermOn_NH_new_mvV(double *qacc,double *qvel,CLT *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double gzi,double *acc_Term,double *acc_Term2,double *vel_Term){
  ABI *abi;

  abi=(ABI *)gcemalloc(sizeof(ABI)*numclt);

  ABAp_prepass_TermOn(clt,qvel,numclt,numatom,crd,vel_Term);

  ABAs_forc(clt,frc,numclt);

  ABAm_mainpass_TermOn_NH_new_mvV(abi,clt,Q,numclt,acc_Term,acc_Term2,vel_Term);

  ABAm_backpass_TermOn_NH_new_mvV(qacc,abi,clt,numclt,qvel,gzi,acc_Term,acc_Term2,vel_Term);

}

double solverABA_NH_chain(double *qacc,double *qvel,CLT *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *zeta_vel,double *zeta_acc,int M,int N,double KBT,double *Q_NH,double tau2,double Temp,double TempB){
  ABI *abi;

  abi=(ABI *)gcemalloc(sizeof(ABI)*numclt);

  ABAp_prepass(clt,qvel,numclt,numatom,crd);
  ABAs_forc(clt,frc,numclt);
  ABAm_mainpass(abi,clt,Q,numclt);
  ABAm_backpass_NH_chain(qacc,abi,clt,numclt,qvel,zeta_vel[0]);
  ABANH_chain_solve(zeta_vel,zeta_acc,Q_NH,M,N,KBT,tau2,Temp,TempB);
}

double solverABA_TermOn_NH_chain(double *qacc,double *qvel,CLT *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *acc_Term,double *acc_Term2,double *vel_Term,double *zeta_vel,double *zeta_acc,int M,int N,double KBT,double *Q_NH,double tau2,double Temp,double TempB){
  ABI *abi;

  abi=(ABI *)gcemalloc(sizeof(ABI)*numclt);

  ABAp_prepass_TermOn(clt,qvel,numclt,numatom,crd,vel_Term);
  ABAs_forc(clt,frc,numclt);
  ABAm_mainpass_TermOn(abi,clt,Q,numclt,acc_Term,acc_Term2,vel_Term);
  ABAm_backpass_TermOn_NH_chain(qacc,abi,clt,numclt,qvel,zeta_vel[0],acc_Term,acc_Term2,vel_Term);
  ABANH_chain_solve(zeta_vel,zeta_acc,Q_NH,M,N,KBT,tau2,Temp,TempB);
}


void ABA_out_formated(FILE *outputfile,double pot,double KE,int i,double dt) {
  fprintf(outputfile,"/***********************************************/\n");
  fprintf(outputfile,"steps            = %d  \n",i);
  fprintf(outputfile,"total time       = %10.3lf ps  \n",dt*(double)i);
  fprintf(outputfile,"toal_energy      = %e kcal/mol  \n",pot+KE);
  fprintf(outputfile,"toal_vertial_energy      = %e kcal/mol  \n",0.0);
  fprintf(outputfile,"\n");
}
