#ifndef INCLUDE_ABAb_INVM
#define INCLUDE_ABAb_INVM

void ABAb_Inverse_mainpass(CLTb *clt,double *qacc,double *qvel,int numclt,int numatom,double *crd);
void ABAb_Inverse_Spacc(double Spacc[6],double qacc,double TMat[6][6],double Coacc[6],double SpaccPart[6]);

#endif








