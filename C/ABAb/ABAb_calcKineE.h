#ifndef INCLUDE_ABAb_KE
#define INCLUDE_ABAb_KE

double ABAb_calcKineE(double *KE,double *KEv,double *PEv,double KEobj,CLTb *clt,double *crd,double *qvel,double q_NVT,double qvel_NVT,double s_NVT,int numclut,int numatom,int MODE);

double ABAb_calcKineE_TermOn(double *KE,double *KEv,double *PEv,double KEobj,CLTb *clt,double *crd,double *qvel,double q_NVT,double qvel_NVT,double s_NVT,double *velTerm,int numclut,int numatom,int MODE,int dof);

double ABAb_calcKineE_new(double *KE,CLTb *clt,double *crd,double *qvel,int numclut,int numatom);

double ABAb_calcKineE_TermOn_new(double *KE,CLTb *clt,double *crd,double *qvel,double q_NVT,double qvel_NVT,double *velTerm,int numclut,int numatom);

double ABAb_calcKineE_TermOn_new_simp(double *KE,CLTb *clt,double *crd,double *qvel,double *velTerm,int numclut,int numatom);

double ABAb_Generate_inivelo(CLTb *clt,double *crd,double *qvel,double *vel_Term,int numclut,int numatom,int TERMMODE, double KbT);

#endif

