#ifndef INCLUDE_ABA_SI
#define INCLUDE_ABA_SI

void ABAs_inertia_matrix(CLT *clt,int nNumClut_all,int num_atom_all,double *crd,double *mass);
void ABAs_Inertia_clust(double Inertia_clust[3][3],
			int num_atom_clust,double *crd_local,double *mass);
void ABAs_InertiaMatrix(double InertiaMatrix[6][6],
			double Inertia_clust[3][3],
			double qcom[3], double *summass,
			int num_atom_clust,
			double *crd_local,double *mass);
int ABA_setJoin(CLT *clt,int nNumClut, int joinflag);

#endif
