#ifndef INCLUDE_ABA_ST
#define INCLUDE_ABA_ST

void ABAs_trans_Matrix(CLT *clt,int nNumClut_all,int num_atom_all,double *crd);

void sub_set_trans_Matrix(double TransMatrix[6][6],
			  int nNumClt,          double trans_A_to_CN[3][3],          int origin_atom_a,
			  int nNumCltminousone, double trans_A_to_CN_minousone[3][3],int origin_atom_a_minousone,
			  double *crd);

#endif
