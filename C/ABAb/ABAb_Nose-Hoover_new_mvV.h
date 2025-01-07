#ifndef INCLUDE_ABAb_NH_new_mvV
#define INCLUDE_ABAb_NH_new_mvV

double GearsConstant5[5];
double Telar_Matrix_NH[5][5];

double ABAbNH_calcKE_new_mvV(double zeta,double s,double s_vel,double Q,double KEobj,double *PEv,double *KEv);
void ABAbNH_set_new_mvV(double s,double s_vel,double zeta,double tau, double *tau2,double *Q, double KEobj,double dt);

void ABAbm_backpass_TermOn_NH_new_mvV(double *qacc,ABIb* abi,CLTb *clt,int numclut,double *qvel,double zeta,double *acc_Term,double *acc_Term2,double *vel_Term);

void ABAbm_backpass_NH_new_mvV(double *qacc,ABIb* abi,CLTb *clt,int numclut,double *qvel,double zeta);
void ABAbm_corBF_NH_TERM_new_mvV(double corABIb[6][6],double corBF[6],double *acc_Term,double *acc_Term2,double *vel_Term);
void ABAbm_mainpass_TermOn_NH_new_mvV(ABIb* abi,CLTb *clt,double *Q,int numclut,double *acc_Term,double *acc_Term2,double *vel_Term);

#endif
