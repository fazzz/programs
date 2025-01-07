#ifndef INCLUDE_ABA_S
#define INCLUDE_ABA_S

void ABAs_local_reference(CLT *clt,int nNumClut_all,int num_atom_all,double *crd);
void sub_trans_A_to_CN(double Mtrans_A_to_CN[3][3], double *crd_local,
		       int nNumCltTar,   int origin_atom_a,
		       int nNumCltCoo,   int terminal_atom_a,
		       int xy_set_atom_a,
		       int nNumAtom,     int num_atom,
		       double *crd);

#endif
