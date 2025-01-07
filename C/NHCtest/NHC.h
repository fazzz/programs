#ifndef INCLUDE_NHC
#define INCLUDE_NHC

double GearsConstant[6];
double Telar_Matrix[6][6];

void NHC_update_pret(double *zeta,double *zeta_vel, double **predict_zeta, double **correct_zeta, int M,double dt);
void NHC_update_cort(double *zeta,double *zeta_vel,double *zeta_acc,double **predict_zeta, double **correct_zeta,int M,double dt);
double NHC_calcKE_new(double *zeta,double *zeta_vel,double *Q,int M,int N,double KBT,double *PEv,double *KEv);
void NHC_set_new(double tau, double *tau2,double *Q, int M,int N,double KBT,double **correct_zeta);
void NHC_solve(double *zeta_vel,double *zeta_acc,double *Q,int M,int N,double KBT,double tau2,double Temp,double TempB);

#endif
