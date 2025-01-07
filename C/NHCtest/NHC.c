
#include <stdio.h>
#include <math.h>

#include "NHC.h"

void NHC_update_pret(double *zeta,double *zeta_vel, double **predict_zeta, double **correct_zeta, int M,double dt) {
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

void NHC_update_cort(double *zeta,double *zeta_vel,double *zeta_acc,double **predict_zeta, double **correct_zeta,int M,double dt) {
  int i,j;

  for (i=0;i<M;++i) {
    for (j=0;j<6;++j) 
      correct_zeta[i][j] = predict_zeta[i][j]+GearsConstant[j]*(0.5*dt*dt*zeta_acc[i]-predict_zeta[i][2]);

    zeta[i]=correct_zeta[i][0];
    zeta_vel[i]=correct_zeta[i][1]/dt;
  }
}

double NHC_calcKE(double *zeta,double *zeta_vel,double *Q,int M,int N,double KBT,double *PEv,double *KEv){
  int i;

  *KEv=0.0;
  for (i=0;i<M;++i) *KEv += 0.5*Q[i]*zeta_vel[i]*zeta_vel[i];

  *PEv = N*KBT*zeta[0];
  for (i=1;i<M;++i) *PEv+=KBT*zeta[i];

}

void NHC_set(double tau, double *tau2,double *Q, int M,int N,double KBT,double **correct_zeta) {
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

  GearsConstant[0] = 3.0/16.0;
  GearsConstant[1] = 251.0/360.0;
  GearsConstant[2] = 1.0;
  GearsConstant[3] = 11.0/18.0;
  GearsConstant[4] = 1.0/6.0;
  GearsConstant[5] = 1.0/60.0;

  for (i=0;i<6;++i){
    for (j=0;j<6;++j){
      if (i != j) Telar_Matrix[i][j] = 0.0;
      else Telar_Matrix[i][i] = 1.0;
    }
  }
	
  Telar_Matrix[0][1] = 1.0;
  Telar_Matrix[0][2] = 1.0;
  Telar_Matrix[0][3] = 1.0;
  Telar_Matrix[0][4] = 1.0;
  Telar_Matrix[0][5] = 1.0;
  Telar_Matrix[1][2] = 2.0;
  Telar_Matrix[1][3] = 3.0;
  Telar_Matrix[1][4] = 4.0;
  Telar_Matrix[1][5] = 5.0;
  Telar_Matrix[2][3] = 3.0;
  Telar_Matrix[2][4] = 6.0;
  Telar_Matrix[2][5] = 10.0;
  Telar_Matrix[3][4] = 4.0;
  Telar_Matrix[3][5] = 10.0;
  Telar_Matrix[4][5] = 5.0;
}

void NHC_solve(double *zeta_vel,double *zeta_acc,double *Q,int M,int N,double KBT,double tau2,double Temp,double TempB) {
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

