#ifndef INCLUDE_ABA_caf
#define INCLUDE_ABA_caf

double ABA_calcAttTorque_H(double *Q, CLT *clt,double *crd,double *q,double *qgoal,double kc,int numclut,FILE *outputfile2,int intervalflag);
double ABA_calcAttTorque_Hwd(double *Q,CLT *clt,double *crd,double *q,double *qvel,double *qgoal,double kc,double kd,int numclut,FILE *outputfile2);
double ABA_calcAttTorque_C(double *Q,CLT *clt,double *crd,double *q,double *qgoal,double kc,int numclut,FILE *outputfile2,double GoalReachedThershold);
void ABA_set_ini(double *q,CLT *clt,double *crd,int numclut);

#endif



