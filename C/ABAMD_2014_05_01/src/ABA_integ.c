#include <stdio.h>
#include <math.h>

//#include "ABA.h" // 2014-06-18
#include "ABAb.h"  // 2014-06-18

void ABA_integ_pret(double *qrot,double *qvel,double *q,double *predict,double *correct,double dt,int numclut) {
  int i,j,k;

  for (i=0;i<numclut;++i) {  
    for (j=0;j<6;++j) predict[i*6+j] = 0.0;
    for (j=0;j<6;++j) 
      for (k=0;k<6;++k) 
	predict[i*6+j] += Telar_Matrix[j][k]*correct[i*6+k];

    qrot[i] = 0.0;  
    for (j=1;j<6;++j) qrot[i] += Telar_Matrix[0][j]*correct[i*6+j];
    qvel[i] = predict[i*6+1]/dt;
  }
}

void ABA_integ_cort(double *qrot,double *qvel,double *q,double *qacc,double *predict,double *correct,double dt,int numclut) {
  int i,j;
  double dqacc;

  for (i=0;i<numclut;++i) {
    dqacc = 0.5*dt*dt*qacc[i]-predict[i*6+2];
  
    for (j=0;j<6;++j) correct[i*6+j] = predict[i*6+j]+GearsConstant[j]*dqacc;
  
    qrot[i] = GearsConstant[0]*dqacc;
    qvel[i] = correct[i*6+1]/dt;
  }
}

void ABA_integ_pret_NVT(double *qvel_NVT,double *q_NVT,double *predict_NVT,double *correct_NVT,double dt) {
  int i,j;

  for (i=0;i<6;++i) predict_NVT[i] = 0.0;
  for (i=0;i<6;++i) 
    for (j=0;j<6;++j) 
      predict_NVT[i] += Telar_Matrix[i][j]*correct_NVT[j];

  *q_NVT = predict_NVT[0];
  *qvel_NVT = predict_NVT[1]/dt;
}

void ABA_integ_cort_NVT(double *qvel_NVT,double *q_NVT,double qacc_NVT,double *predict_NVT,double *correct_NVT,double dt) {
  int i,j;
  double dqacc;

  dqacc = 0.5*dt*dt*qacc_NVT-predict_NVT[2];
  for (i=0;i<6;++i) correct_NVT[i] = predict_NVT[i]+GearsConstant[i]*dqacc;  
  *q_NVT = correct_NVT[0];
  *qvel_NVT = correct_NVT[1]/dt;
}

void ABA_integ_set(double *q,double *qvel,double *predict,double *correct,int numclut,double dt) {
  int i,j;

  for (i=0;i<numclut;++i){
    for (j=0;j<6;++j){
      correct[i*6+j]=0.0;
      predict[i*6+j]=0.0;
    }
  }
  for (i=0;i<numclut;++i){
    correct[i*6+0]=q[i];
    correct[i*6+1]=qvel[i]*dt;
  }

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

void ABA_integ_set_NVT(double q_NVT,double qvel_NVT,double *predict_NVT,double *correct_NVT,double dt) {
  int i;

  for (i=0;i<6;++i){
    correct_NVT[i]=0.0;
    predict_NVT[i]=0.0;
  }
  correct_NVT[0]=q_NVT;
  predict_NVT[0]=q_NVT;
  correct_NVT[1]=qvel_NVT*dt;
  predict_NVT[2]=qvel_NVT*dt;
}

void ABA_integ_pret_Term(double **predict_Term,double **predict_Term3,double **correct_Term,double **correct_Term3,double *vel_Term,double *delta_Term,double dt){
  int i,j,k;

  for (i=0;i<6;++i) delta_Term[i] = 0.0;

  for (i=0;i<6;++i) {
    for (j=0;j<6;++j) {
      predict_Term[i][j] = 0.0;
      predict_Term3[i][j] = 0.0;
    }
  }
  for (i=0;i<6;++i) {
    for (j=0;j<6;++j) {
      for (k=0;k<6;++k) {
  	predict_Term[i][j] += Telar_Matrix[j][k]*correct_Term[i][k];
  	predict_Term3[i][j] += Telar_Matrix[j][k]*correct_Term3[i][k];
      }
    }
  }
  for (i=0;i<6;++i) {
    vel_Term[i]=predict_Term[i][1]/dt;
    for (j=1;j<6;++j) delta_Term[i] += Telar_Matrix[i][j]*correct_Term[i][j];
  }

}

void ABA_integ_cort_Term(double **predict_Term,double **predict_Term3,double **correct_Term,double **correct_Term3,double *acc_Term,double *acc_Term2,double *vel_Term,double *delta_Term,double dt) {
  int i,j,k;
  double delta_acc_Term[6],delta_acc_Term3[6];
  
  for (i=0;i<6;++i) {
    delta_acc_Term[i] = 0.5*dt*dt*acc_Term[i]-predict_Term[i][2];
    delta_acc_Term3[i] = 0.5*dt*dt*acc_Term2[i]-predict_Term3[i][2];
  }
  for (i=0;i<6;++i) {
    for (j=0;j<6;++j) {
      correct_Term[i][j] = 0.0;
      correct_Term3[i][j] = 0.0;
    }
  }
  for (i=0;i<6;++i) {
    for (j=0;j<6;++j) {
      correct_Term[i][j] = predict_Term[i][j]+GearsConstant[j]*delta_acc_Term[i];
      correct_Term3[i][j] = predict_Term3[i][j]+GearsConstant[j]*delta_acc_Term3[i];
    }
  }
  
  for (i=0;i<6;++i) {
    vel_Term[i]=correct_Term[i][1]/dt;
    delta_Term[i]=GearsConstant[0]*delta_acc_Term[i];
  }
}

