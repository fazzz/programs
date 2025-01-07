#ifndef INCLUDE_ABAb_ST
#define INCLUDE_ABAb_ST

void ABAbs_trans_Matrix(CLTb *clt,int nNumClut_all,int num_atom_all,double *crd);

void sub_set_trans_Matrix(double TransMatrix[6][6],
			  int nNumClt,          double trans_A_to_CN[3][3],          int origin_atom_a,
			  int nNumCltminousone, double trans_A_to_CN_minousone[3][3],int origin_atom_a_minousone,
			  double *crd);

#endif
