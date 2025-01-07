#ifndef INCLUDE_ABAb_SI
#define INCLUDE_ABAb_SI

void ABAbs_inertia_matrix(CLTb *clt,int nNumClut_all,int num_atom_all,double *crd,double *mass);
void ABAbs_Inertia_clust(double Inertia_clust[3][3],
			int num_atom_clust,double *crd_local,double *mass);
void ABAbs_InertiaMatrix(double InertiaMatrix[6][6],
			double Inertia_clust[3][3],
			double qcom[3], double *summass,
			int num_atom_clust,
			double *crd_local,double *mass);
int ABAb_setJoin(CLTb *clt,int nNumClut, int joinflag);

#endif
