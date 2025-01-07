#ifndef INCLUDE_ABA_NH
#define INCLUDE_ABA_NH

double GearsConstant_NH[5];
double Telar_Matrix_NH[5][5];

void ABANH_update_pret(double *s,double *s_vel,double predict_s[6], double correct_s[6],double dt);
void ABANH_update_cort(double *gzi, double *s, double *s_vel, double s_acc, double predict_s[6], double correct_s[6], double dt);

void ABANH_calD(double s, double *s_acc, double *s_vel,double Temp, double TempB, double tau2);
double ABANH_calcKE(double s,double s_vel,double Q,double KEobj,double *PEv,double *KEv);

void ABANH_set(double s,double s_vel,double gzi,double predict_s[6],double correct_s[6],double tau, double *tau2, double *Q, double KEobj,double dt);

void ABAm_backpass_TermOn_NH(double *qacc,ABI* abi,CLT *clt,int numclut,double *qvel,double s,double s_vel,double *s_acc,double tau2,double *acc_Term,double *acc_Term2,double *vel_Term,double Temp,double TempB);

void ABAm_backpass_NH(double *qacc,ABI* abi,CLT *clt,int numclut,double *qvel,double s,double s_vel,double *s_acc,double tau2,double Temp,double TempB);

void ABAm_corBF_NH_TERM(double corABI[6][6],double corBF[6],double *acc_Term,double *acc_Term2,double *vel_Term,double s,double s_vel);

void ABAm_mainpass_TermOn_NH(ABI* abi,CLT *clt,double *Q,int numclut,double *acc_Term,double *acc_Term2,double *vel_Term,double s,double s_vel);

#endif
