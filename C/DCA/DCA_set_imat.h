#ifndef INCLUDE_DCA_SI
#define INCLUDE_DCA_SI

void DCAs_inertia_matrix(CLT *clt,int nNumClut_all,int num_atom_all,double *crd,double *mass);
void DCAs_Inertia_clust(double Inertia_clust[3][3],
			int num_atom_clust,double *crd_local,double *mass);
void DCAs_InertiaMatrix(double InertiaMatrix[6][6],
			double Inertia_clust[3][3],
			int num_atom_clust,
			double *crd_local,double *mass);

#endif
