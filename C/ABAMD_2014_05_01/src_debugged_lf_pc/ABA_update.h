#ifndef INCLUDE_ABA_U
#define INCLUDE_ABA_U

void ABA_update(CLT *clt,double *crd,double *deltaq,int numclut,int numatom);
void ABA_update_quaternion(double delta_dihed,double *crd, CLT *clt,
			   int nNumClt,int nNumAtomALL, int nNumClutLast, int nNumCltPt);
//void ABA_update_Term(double *crd,double *delta, int numatom, CLT *clt); // 2014-06-19
//void ABA_update_Term(double *crd,double *delta, int numatom, CLT *clt, int numclut); // 2014-06-24
//void ABA_update_Term(double *crd,double *delta, int numatom, CLT *clt, int numclut, double trans_A_to_CN_terminal[3][3],double l_Term[3]);
//void ABA_update_Term(double *crd,double *delta, int numatom, CLT *clt, int numclut,   // 2014-06-30 // 2014-08-13
//		     double *trans_A_to_CN_terminal/*[3][3]*/,double *l_Term/*[3]*/); // 2014-06-30 // 2014-08-13
void ABA_update_Term(double *crd,     // 2014-08-13
		     double *delta, // 2014-08-13
		     //		     double delta[6], // 2014-08-13
		     int numatom, CLT *clt, int numclut,  // 2014-08-13
		     double *trans_A_to_CN_terminal/*[3][3]*/,double *l_Term/*[3]*/); // 2014-08-13
void ABA_update_Term2(double *crd,double *crd_Term,int numatom,              // 2014-06-30
		      double *trans_A_to_CN_terminal,double *l_Term);                 // 2014-06-30

#endif

