#ifndef INCLUDE_ABAb_I
#define INCLUDE_ABAb_I

double GearsConstant[6];
double Telar_Matrix[6][6];

void ABAb_integ_set(double *q,double *qvel,double *predict,double *correct,int numclut,double dt);

void ABAb_integ_pret(double *qrot,double *qvel,double *q,double *predict,double *correct,double dt,int numclut);
void ABAb_integ_cort(double *qrot,double *qvel,double *q,double *qacc,double *predict,double *correct,double dt,int numclut);


void ABAb_integ_set_NVT(double q_NVT,double qvel_NVT,double *predict_NVT,double *correct_NVT,double dt);

void ABAb_integ_pret_NVT(double *qvel_NVT,double *q_NVT,double *predict_NVT,double *correct_NVT,double dt);
void ABAb_integ_cort_NVT(double *qvel_NVT,double *q_NVT,double qacc_NVT,double *predict_NVT,double *correct_NVT,double dt);

void ABAb_integ_pret_Term(double **predict_Term,double **predict_Term3,double **correct_Term,double **correct_Term3,double *vel_Term,double *delta_Term,double dt);
void ABAb_integ_cort_Term(double **predict_Term,double **predict_Term3,double **correct_Term,double **correct_Term3,double *acc_Term,double *acc_Term2,double *vel_Term,double *delta_Term,double dt);
  
#endif

