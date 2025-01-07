
#ifndef INCLUDE_TACCM
#define INCLUDE_TACCM

double GearsConstantZ[6],GearsConstantZ5[5];
double Telar_MatrixZ[6][6];

double TACCM_calc_eff_FF(double *theta, double *Z,  int numZ,double Kapa, double *Q, int **pairs,double pi);

double TACCM_calc_eff_FF_Z(double *Z,int numZ,double *theta,double KZ,double *f, double pi);

double TACCM_calc_eff_FF_Z_2(double *Z,int numZ,double *theta,double KZ,double *f, double *PE,double pi);

//double TACCM_MD_Propagetor_NH_MP1998_Z(double *Z,double *velZ,double massZ,double *theta,
//				       double *zeta,double *V_zeta,double Q,double NfKT,int numZ,
//				       double *KEv,double *PEv,
//				       double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
//				       double KZ,double *PEZ,double *f,double pi);


//double TACCM_MD_Propagetor_NH_Single_part_MP1996(double *vel,double massZ,double *zeta,
//						 double *V_zeta,double Q,double NfKT,
//						 int numZ,int nc,double wdt4[3],double wdt2[3]);

//double TACCM_MD_Propagetor_vV_NVE_Z(double *Z, double *velZ,double massZ,double *theta,int numZ,
//				    double dt, double KZ,double *PEZ,double *f,double pi);

double TACCM_CTheta(double *crd,int numatom,double *theta, int numdihe, int **pairs, double pi);

double TACCM_MD_Generate_inivelo(double *velZ,double mass_Z,int numZ,double KbT);

void TACCM_integ_pret_Z(double **predict_Z,double **correct_Z,double *Z,double *vel_Z,int numZ,double dt,double pi);

void TACCM_integ_cort_Z(double **predict_Z,double **correct_Z,double *acc_Z,double *Z,double *vel_Z,int numZ,double dt,double pi);

double TACCM_calcKineE_Z(double *KE,double massZ,double *vel_Z,int numZ);

double TACCM_solver_NH_Z(double *accZ,double *velZ,double massZ,double *frcZ,int numZ,double gzi,double *gzi_vel,double tau2,double Temp,double TempB);

void TACCM_NH_update_pret_new(double *gzi,double *gzi_vel,double predict_gzi[5], double correct_gzi[5],
			      double *s,double *s_vel,double predict_s[5], double correct_s[5],double dt);

void TACCM_NH_update_cort_new(double *gzi, double gzi_vel, double *s, double *s_vel,
			      double predict_gzi[5], double correct_gzi[5],
			      double predict_s[5], double correct_s[5], double dt);

double TACCM_NH_calcKE_new(double gzi,double s,double s_vel,double Q,double KEobj,double *PEv,double *KEv);

void TACCM_NH_set_new(double s,double s_vel,double gzi,double predict_gzi[5],double correct_gzi[5],
		      double predict_s[5],double correct_s[5],double tau, double *tau2, 
		      double *Q, double KEobj,double dt);

#endif
