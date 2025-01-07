#ifndef INCLUDE_ABA_NH_new_mvV
#define INCLUDE_ABA_NH_new_mvV

double GearsConstant5[5];
double Telar_Matrix_NH[5][5];

double ABANH_calcKE_new_mvV(double zeta,double s,double s_vel,double Q,double KEobj,double *PEv,double *KEv);
void ABANH_set_new_mvV(double s,double s_vel,double zeta,double tau, double *tau2,double *Q, double KEobj,double dt);

void ABAm_backpass_TermOn_NH_new_mvV(double *qacc,ABI* abi,CLT *clt,int numclut,double *qvel,double zeta,double *acc_Term,double *acc_Term2,double *vel_Term);

void ABAm_backpass_NH_new_mvV(double *qacc,ABI* abi,CLT *clt,int numclut,double *qvel,double zeta);
void ABAm_corBF_NH_TERM_new_mvV(double corABI[6][6],double corBF[6],double *acc_Term,double *acc_Term2,double *vel_Term);
void ABAm_mainpass_TermOn_NH_new_mvV(ABI* abi,CLT *clt,double *Q,int numclut,double *acc_Term,double *acc_Term2,double *vel_Term);

#endif
