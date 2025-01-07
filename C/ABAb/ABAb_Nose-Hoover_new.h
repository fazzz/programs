#ifndef INCLUDE_ABAb_NH_new
#define INCLUDE_ABAb_NH_new

double GearsConstant5[5];
double Telar_Matrix_NH[5][5];

void ABAbNH_update_pret_new(double *gzi,double predict_gzi[5], double correct_gzi[5],
			   double *s,double *s_vel,double predict_s[5], double correct_s[5],double dt);
void ABAbNH_update_cort_new(double *gzi, double gzi_vel, double *s, double *s_vel,
			   double predict_gzi[5], double correct_gzi[5],
			   double predict_s[5], double correct_s[5], double dt);
double ABAbNH_calcKE_new(double gzi,double s,double s_vel,double Q,double KEobj,double *PEv,double *KEv);
void ABAbNH_set_new(double s,double s_vel,double gzi,double predict_gzi[5],double correct_gzi[5],
		   double predict_s[5],double correct_s[5],double tau, double *tau2, 
		   double *Q, double KEobj,double dt);
void ABAbm_backpass_TermOn_NH_new(double *qacc,ABIb* abi,CLTb *clt,int numclut,double *qvel,double gzi,double *gzi_vel,double tau2,double *acc_Term,double *acc_Term2,double *vel_Term,double Temp,double TempB);
void ABAbm_backpass_NH_new(double *qacc,ABIb* abi,CLTb *clt,int numclut,double *qvel,double gzi,double *gzi_vel,double tau2,double Temp,double TempB);
void ABAbm_corBF_NH_TERM_new(double corABIb[6][6],double corBF[6],double *acc_Term,double *acc_Term2,double *vel_Term);
void ABAbm_mainpass_TermOn_NH_new(ABIb* abi,CLTb *clt,double *Q,int numclut,double *acc_Term,double *acc_Term2,double *vel_Term);

#endif
