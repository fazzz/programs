#ifndef INCLUDE_ABAb_B
#define INCLUDE_ABAb_B

void ABAbm_backpass(double *qacc,ABIb* abi,CLTb *clt,int numclut,double *qvel,double q_NVT,double qvel_NVT,double *qacc_NVT,double s_NVT,double KE,double KEobj,int MODE);
void ABAbm_backpass_TermOn(double *qacc,ABIb* abi,CLTb *clt,int numclut,double *qvel,double q_NVT,double qvel_NVT,double *qacc_NVT,double s_NVT,double *acc_Term,double *acc_Term2,double KE,double KEobj,int MODE);
void ABAbm_backpass_TERM(double acc_Term2[6],double acc_Term[6],double vel_Term[6],double s_NVT,double svel_NVT,int MODE);

#endif
