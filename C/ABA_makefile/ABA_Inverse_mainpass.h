#ifndef INCLUDE_ABA_INVM
#define INCLUDE_ABA_INVM

void ABA_Inverse_mainpass(CLT *clt,double *qacc,double *qvel,int numclt,int numatom,double *crd);
void ABA_Inverse_Spacc(double Spacc[6],double qacc,double TMat[6][6],double Coacc[6],double SpaccPart[6]);

#endif








