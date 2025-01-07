#ifndef INCLUDE_ABAb_caf
#define INCLUDE_ABAb_caf

double ABAb_calcAttTorque_H(double *Q, CLTb *clt,double *crd,double *q,double *qgoal,double kc,int numclut,FILE *outputfile2,int intervalflag);
double ABAb_calcAttTorque_Hwd(double *Q,CLTb *clt,double *crd,double *q,double *qvel,double *qgoal,double kc,double kd,int numclut,FILE *outputfile2);
double ABAb_calcAttTorque_C(double *Q,CLTb *clt,double *crd,double *q,double *qgoal,double kc,int numclut,FILE *outputfile2,double GoalReachedThershold);
void ABAb_set_ini(double *q,CLTb *clt,double *crd,int numclut);

#endif



