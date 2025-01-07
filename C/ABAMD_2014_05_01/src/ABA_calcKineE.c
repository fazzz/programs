#include <stdio.h>
#include <math.h>

//#include "ABA.h" // 2014-06-18
#include "ABAb.h"  // 2014-06-18
#include "LA.h"
#include "EF.h"

#include "RAND.h"
#include "BOXMULL.h"

#define ON 1
#define OFF 0

double ABA_calcKineE(double *KE,double *KEv,double *PEv,double KEobj,CLT *clt,double *crd,double *qvel,double q_NVT,double qvel_NVT,double s_NVT,int numclut,int numatom,int MODE){
  int i,j,k;
  double *vec,*mat;
  double UNIT=418.4070;

  ABAp_prepass(clt,qvel,numclut,numatom,crd);

  *KE=0.0;
  *KEv=0.0;
  *PEv=0.0;
  vec=(double *)gcemalloc(sizeof(double)*6);
  mat=(double *)gcemalloc(sizeof(double)*6*6);

  for (i=1;i<numclut;++i) {
    for (j=0;j<6;++j) {
      vec[j]=clt[i].Spvel[j];
      for (k=0;k<6;++k)
	mat[j*6+k]=clt[i].IM[j][k];
    }
    //    *KE+=0.5*vtmvmult(vec,mat,vec,6)/UNIT;
    *KE+=vtmvmult(vec,mat,vec,6);
  }
  *KE=0.5*(*KE)/(4.18407*100.0)/*UNIT*/;

  if (MODE == NVT) {
    *KEv = 0.5*s_NVT*(qvel_NVT/q_NVT)*(qvel_NVT/q_NVT)/UNIT;
    *PEv = KEobj*log(q_NVT);
  }

  return 0.0; // 2014-07-04
}

double ABA_calcKineE_TermOn(double *KE,double *KEv,double *PEv,double KEobj,CLT *clt,double *crd,double *qvel,double q_NVT,double qvel_NVT,double s_NVT,double *velTerm,int numclut,int numatom,int MODE,int dof){
  int i,j,k;
  double *vec,*mat;
  double UNIT=418.4070;

  ABAp_prepass_TermOn(clt,qvel,numclut,numatom,crd,velTerm);

  *KE=0.0;
  *KEv=0.0;
  *PEv=0.0;
  vec=(double *)gcemalloc(sizeof(double)*6);
  mat=(double *)gcemalloc(sizeof(double)*6*6);

  for (i=0;i<numclut;++i) {
    for (j=0;j<6;++j) {
      vec[j]=clt[i].Spvel[j];
      for (k=0;k<6;++k)
	mat[j*6+k]=clt[i].IM[j][k];
    }
    //    *KE+=0.5*vtmvmult(vec,mat,vec,6)/UNIT;
    *KE+=vtmvmult(vec,mat,vec,6);
  }
  *KE=0.5*(*KE)/(4.18407*100.0)/*/UNIT*/;

  if (MODE == NVT) {
    *KEv = 0.5*s_NVT*(qvel_NVT/q_NVT)*(qvel_NVT/q_NVT)/UNIT;
    *PEv = KEobj*log(q_NVT);
  }

  return 0.0; // 2104-07-04
}

double ABA_calcKineE_new(double *KE,CLT *clt,double *crd,double *qvel,int numclut,int numatom){
  int i,j,k;
  double *vec,*mat;
  double UNIT=418.4070;

  ABAp_prepass(clt,qvel,numclut,numatom,crd);

  *KE=0.0;
  vec=(double *)gcemalloc(sizeof(double)*6);
  mat=(double *)gcemalloc(sizeof(double)*6*6);

  for (i=1;i<numclut;++i) {
    for (j=0;j<6;++j) {
      vec[j]=clt[i].Spvel[j];
      for (k=0;k<6;++k)
	mat[j*6+k]=clt[i].IM[j][k];
    }
    //    *KE+=0.5*vtmvmult(vec,mat,vec,6)/UNIT;
    *KE+=vtmvmult(vec,mat,vec,6);
  }
  *KE=0.5*(*KE)/(4.18407*100.0)/*UNIT*/;

  return 0.0; // 2014-07-04
}

double ABA_calcKineE_TermOn_new(double *KE,CLT *clt,double *crd,double *qvel,double q_NVT,double qvel_NVT,double *velTerm,int numclut,int numatom){
  int i,j,k;
  double *vec,*mat;
  double UNIT=418.4070;

  ABAp_prepass_TermOn(clt,qvel,numclut,numatom,crd,velTerm);

  *KE=0.0;
  vec=(double *)gcemalloc(sizeof(double)*6);
  mat=(double *)gcemalloc(sizeof(double)*6*6);

  for (i=0;i<numclut;++i) {
    for (j=0;j<6;++j) {
      vec[j]=clt[i].Spvel[j];
      for (k=0;k<6;++k)
	mat[j*6+k]=clt[i].IM[j][k];
    }
    //    *KE+=0.5*vtmvmult(vec,mat,vec,6)/UNIT;
    *KE+=vtmvmult(vec,mat,vec,6);
  }
  *KE=0.5*(*KE)/(4.18407*100.0)/*/UNIT*/;

  return 0.0; // 2014-07-04
}

double ABA_calcKineE_TermOn_new_simp(double *KE,CLT *clt,double *crd,double *qvel,double *velTerm,int numclut,int numatom){
  int i,j,k;
  double *vec,*mat;
  double UNIT=418.4070;

  ABAp_prepass_TermOn(clt,qvel,numclut,numatom,crd,velTerm);

  *KE=0.0;
  vec=(double *)gcemalloc(sizeof(double)*6);
  mat=(double *)gcemalloc(sizeof(double)*6*6);

  for (i=0;i<numclut;++i) {
    for (j=0;j<6;++j) {
      vec[j]=clt[i].Spvel[j];
      for (k=0;k<6;++k)
	mat[j*6+k]=clt[i].IM[j][k];
    }
    *KE+=vtmvmult(vec,mat,vec,6);
  }
  *KE=0.5*(*KE)/(4.18407*100.0);

  return 0.0; // 2014-07-04
}


double ABA_Generate_inivelo(CLT *clt,double *crd,double *qvel,double *vel_Term,int numclut,int numatom,int TERMMODE, double KbT) {
  int i,j,k;
  double *vec,*mat;
  double KE=0.0;
  double UNITT=418.4070;
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  if (TERMMODE==ON)
    for (i=0;i<6;++i) 
      vel_Term[i]=Box_Muller(i,0.0,KbT/clt[0].sum_mass);

  for (i=1;i<numclut;++i) 
    qvel[i]=Box_Muller(i,0.0,KbT/clt[i].I[2][2]);

  KE=0.0;
  vec=(double *)gcemalloc(sizeof(double)*6);
  mat=(double *)gcemalloc(sizeof(double)*6*6);

  ABAp_prepass_TermOn(clt,qvel,numclut,numatom,crd,vel_Term);

  for (i=0;i<numclut;++i) {
    for (j=0;j<6;++j) {
      vec[j]=clt[i].Spvel[j];
      for (k=0;k<6;++k)
	mat[j*6+k]=clt[i].IM[j][k];
    }
    KE+=vtmvmult(vec,mat,vec,6);
  }
  KE=0.5*KE/(4.18407*100.0);

  return KE;
}
