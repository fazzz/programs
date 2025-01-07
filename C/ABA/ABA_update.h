#ifndef INCLUDE_ABA_U
#define INCLUDE_ABA_U

void ABA_update(CLT *clt,double *crd,double *deltaq,int numclut,int numatom);
void ABA_update_quaternion(double delta_dihed,double *crd, CLT *clt,
			   int nNumClt,int nNumAtomALL, int nNumClutLast, int nNumCltPt);

#endif

