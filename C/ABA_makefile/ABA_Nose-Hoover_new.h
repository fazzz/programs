#ifndef INCLUDE_ABA_NH_new
#define INCLUDE_ABA_NH_new

double GearsConstant5[5];
double Telar_Matrix_NH[5][5];

void ABANH_update_pret_new(double *gzi,double predict_gzi[5], double correct_gzi[5],
			   double *s,double *s_vel,double predict_s[5], double correct_s[5],double dt);
void ABANH_update_cort_new(double *gzi, double gzi_vel, double *s, double *s_vel,
			   double predict_gzi[5], double correct_gzi[5],
			   double predict_s[5], double correct_s[5], double dt);
double ABANH_calcKE_new(double gzi,double s,double s_vel,double Q,double KEobj,double *PEv,double *KEv);
void ABANH_set_new(double s,double s_vel,double gzi,double predict_gzi[5],double correct_gzi[5],
		   double predict_s[5],double correct_s[5],double tau, double *tau2, 
		   double *Q, double KEobj,double dt);
void ABAm_backpass_TermOn_NH_new(double *qacc,ABI* abi,CLT *clt,int numclut,double *qvel,double gzi,double *gzi_vel,double tau2,double *acc_Term,double *acc_Term2,double *vel_Term,double Temp,double TempB);
void ABAm_backpass_NH_new(double *qacc,ABI* abi,CLT *clt,int numclut,double *qvel,double gzi,double *gzi_vel,double tau2,double Temp,double TempB);
void ABAm_corBF_NH_TERM_new(double corABI[6][6],double corBF[6],double *acc_Term,double *acc_Term2,double *vel_Term);
void ABAm_mainpass_TermOn_NH_new(ABI* abi,CLT *clt,double *Q,int numclut,double *acc_Term,double *acc_Term2,double *vel_Term);

#endif
