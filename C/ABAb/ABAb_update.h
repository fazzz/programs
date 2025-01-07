#ifndef INCLUDE_ABAb_U
#define INCLUDE_ABAb_U

void ABAb_update(CLTb *clt,double *crd,double *deltaq,int numclut,int numatom);
void ABAb_update_quaternion(double delta_dihed,double *crd,CLTb *clt,
			    int nNumClt,int nNumAtomALL, int nNumClutLast, int nNumCltPt,
			    int nNumAtomP, int nNumAtomO);

#endif

