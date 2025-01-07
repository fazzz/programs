
#include <stdio.h>
#include <math.h>

//#include "ABA.h" // 2014-06-18
#include "ABAb.h"  // 2014-06-18

void ABANH_chain_update_pret(double *zeta,double *zeta_vel, double **predict_zeta, double **correct_zeta, int M,double dt) {
  int i,j,k;

  for (i=0;i<M;++i)
    for (j=0;j<6;++j) 
      predict_zeta[i][j] = 0.0;

  for (i=0;i<M;++i) {
    for (j=0;j<6;++j) {
      for (k=0;k<6;++k) {
	predict_zeta[i][j] += Telar_Matrix[j][k]*correct_zeta[i][k];
      }
    }

    zeta[i]=predict_zeta[i][0];
    zeta_vel[i]=predict_zeta[i][1]/dt;
  }
}

void ABANH_chain_update_cort(double *zeta,double *zeta_vel,double *zeta_acc,double **predict_zeta, double **correct_zeta,int M,double dt) {
  int i,j;

  for (i=0;i<M;++i) {
    for (j=0;j<6;++j) 
      correct_zeta[i][j] = predict_zeta[i][j]+GearsConstant[j]*(0.5*dt*dt*zeta_acc[i]-predict_zeta[i][2]);

    zeta[i]=correct_zeta[i][0];
    zeta_vel[i]=correct_zeta[i][1]/dt;
  }
}

double ABANH_chain_calcKE_new(double *zeta,double *zeta_vel,double *Q,int M,int N,double KBT,double *PEv,double *KEv){
  int i;

  *KEv=0.0;
  for (i=0;i<M;++i) *KEv += 0.5*Q[i]*zeta_vel[i]*zeta_vel[i];

  *PEv = N*KBT*zeta[0];
  for (i=1;i<M;++i) *PEv+=KBT*zeta[i];

  return 0.0; // 2014-07-04
}

void ABANH_chain_set_new(double tau, double *tau2,double *Q, int M,int N,double KBT,double **correct_zeta) {
  int i,j;
  double UNIT=418.4070;
  double pi;

  pi=acos(-1.0);

  tau=tau/2.0/pi;
  *tau2=tau*tau; 
  Q[0]=(*tau2)*N*KBT;   

  for (i=1;i<M;++i) Q[i]=(*tau2)*KBT;   

  for (i=0;i<M;++i) 
    for (j=0;j<6;++j) 
      correct_zeta[i][j]=0.0;

}

void ABANH_chain_solve(double *zeta_vel,double *zeta_acc,double *Q,int M,int N,double KBT,double tau2,double Temp,double TempB) {
  int i,j,k;
  double UNIT=418.4070;

  zeta_acc[0] = 1.0/(tau2)*(Temp/TempB-1.0);

  if (M>1) {
    zeta_acc[0]-=zeta_vel[0]*zeta_vel[1];
    for (i=1;i<M-1;++i) {
      zeta_acc[i] = 1.0/Q[i]*(Q[i-1]*zeta_vel[i-1]*zeta_vel[i-1]-KBT)-zeta_vel[i]*zeta_vel[i+1];
    }
    zeta_acc[i] = 1.0/Q[i]*(Q[i-1]*zeta_vel[i-1]*zeta_vel[i-1]-KBT);
  }

}

void ABAm_backpass_TermOn_NH_chain(double *qacc,ABI* abi,CLT *clt,int numclut,double *qvel,double zeta_vel,double *acc_Term,double *acc_Term2,double *vel_Term) {
  int i,j,k;
  double add[6];
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

  add[0]=0.0;
  add[1]=0.0;
  add[2]=0.0;
  add[3]=(-vel_Term[2]*vel_Term[3+1]+vel_Term[1]*vel_Term[3+2]);
  add[4]=( vel_Term[2]*vel_Term[3+0]-vel_Term[0]*vel_Term[3+2]);
  add[5]=(-vel_Term[1]*vel_Term[3+0]+vel_Term[0]*vel_Term[3+1]);

  for(i=0;i<6;++i) { 
    acc_Term2[i] -= zeta_vel*vel_Term[i];
    acc_Term[i]=acc_Term2[i]-add[i];
  }

  for(i=0;i<numclut;++i) qacc[i] -= zeta_vel*qvel[i];
}

void ABAm_backpass_NH_chain(double *qacc,ABI* abi,CLT *clt,int numclut,double *qvel,double zeta_vel) {
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

  for(i=0;i<numclut;++i) qacc[i] -= zeta_vel*qvel[i];
}


