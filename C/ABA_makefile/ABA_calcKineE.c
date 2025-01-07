#include <stdio.h>
#include <math.h>

#include "ABA.h"
#include "LA.h"
#include "EF.h"

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

  for (i=0;i<numclut;++i) {
    for (j=0;j<6;++j) {
      vec[j]=clt[i].Spvel[j];
      for (k=0;k<6;++k)
	mat[j*6+k]=clt[i].IM[j][k];
    }
    *KE+=0.5*vtmvmult(vec,mat,vec,6)/UNIT;
  }

  if (MODE == NVT) {
    *KEv = 0.5*s_NVT*(qvel_NVT/q_NVT)*(qvel_NVT/q_NVT)/UNIT;
    *PEv = KEobj*log(q_NVT);
  }
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
    *KE+=0.5*vtmvmult(vec,mat,vec,6)/UNIT;
  }

  if (MODE == NVT) {
    *KEv = 0.5*s_NVT*(qvel_NVT/q_NVT)*(qvel_NVT/q_NVT)/UNIT;
    *PEv = KEobj*log(q_NVT);
  }
}

