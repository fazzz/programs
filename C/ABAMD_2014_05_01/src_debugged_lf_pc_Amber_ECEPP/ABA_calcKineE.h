#ifndef INCLUDE_ABA_KE
#define INCLUDE_ABA_KE

double ABA_calcKineE(double *KE,double *KEv,double *PEv,double KEobj,CLT *clt,double *crd,double *qvel,double q_NVT,double qvel_NVT,double s_NVT,int numclut,int numatom,int MODE);

double ABA_calcKineE_TermOn(double *KE,double *KEv,double *PEv,double KEobj,CLT *clt,double *crd,double *qvel,double q_NVT,double qvel_NVT,double s_NVT,double *velTerm,int numclut,int numatom,int MODE,int dof);

double ABA_calcKineE_new(double *KE,CLT *clt,double *crd,double *qvel,int numclut,int numatom);

double ABA_calcKineE_TermOn_new(double *KE,CLT *clt,double *crd,double *qvel,double q_NVT,double qvel_NVT,double *velTerm,int numclut,int numatom);

double ABA_calcKineE_TermOn_new_simp(double *KE,CLT *clt,double *crd,double *qvel,double *velTerm,int numclut,int numatom);

double ABA_Generate_inivelo(CLT *clt,double *crd,double *qvel,double *vel_Term,int numclut,int numatom,int TERMMODE, double KbT);

#endif

