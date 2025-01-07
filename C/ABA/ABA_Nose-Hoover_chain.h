#ifndef INCLUDE_ABA_NH_chain
#define INCLUDE_ABA_NH_chain

void ABANH_chain_update_pret(double *zeta,double *zeta_vel, double **predict_zeta, double **correct_zeta, int M,double dt);

void ABANH_chain_update_cort(double *zeta,double *zeta_vel,double *zeta_acc,double **predict_zeta, double **correct_zeta,int M,double dt);

double ABANH_chain_calcKE_new(double *zeta,double *zeta_vel,double *Q,int M,int N,double KBT,double *PEv,double *KEv);

void ABANH_chain_set_new(double tau, double *tau2,double *Q, int M,int N,double KBT,double **correct_zeta);

void ABANH_chain_solve(double *zeta_vel,double *zeta_acc,double *Q,int M,int N,double KBT,double tau2,double Temp,double TempB);

void ABAm_backpass_TermOn_NH_chain(double *qacc,ABI* abi,CLT *clt,int numclut,double *qvel,double zeta_vel,double *acc_Term,double *acc_Term2,double *vel_Term);

void ABAm_backpass_NH_chain(double *qacc,ABI* abi,CLT *clt,int numclut,double *qvel,double zeta_vel);

#endif
