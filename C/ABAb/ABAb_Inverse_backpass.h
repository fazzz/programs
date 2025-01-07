#ifndef INCLUDE_ABAb_INVB
#define INCLUDE_ABAb_INVB

void ABAb_Inverse_backpass(ABIb* abi,CLTb *clt,double *Q,int numclut);
double ABAbb_cQ(double F[6], double IM[6][6], double Spacc[6],double Cofrc[6],double Spfrc[6], double TM[4][6][6],double frc[4][6],int num_branch);

#endif



