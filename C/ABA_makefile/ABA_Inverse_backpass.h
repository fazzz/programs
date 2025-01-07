#ifndef INCLUDE_ABA_INVB
#define INCLUDE_ABA_INVB

void ABA_Inverse_backpass(ABI* abi,CLT *clt,double *Q,int numclut);
double ABAb_cQ(double F[6], double IM[6][6], double Spacc[6],double Cofrc[6],double Spfrc[6], double TM[4][6][6],double frc[4][6],int num_branch);

#endif



