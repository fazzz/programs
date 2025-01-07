#ifndef INCLUDE_ABA_KE
#define INCLUDE_ABA_KE

double ABA_calcKineE(double *KE,double *KEv,double *PEv,double KEobj,CLT *clt,double *crd,double *qvel,double q_NVT,double qvel_NVT,double s_NVT,int numclut,int numatom,int MODE);

double ABA_calcKineE_TermOn(double *KE,double *KEv,double *PEv,double KEobj,CLT *clt,double *crd,double *qvel,double q_NVT,double qvel_NVT,double s_NVT,double *velTerm,int numclut,int numatom,int MODE,int dof);

#endif

