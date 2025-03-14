

//int calc_dPC(int numdihed, int numstep, int numpc);
//int scandtraj(FILE *inputfile, int numdihed, int numstep);

int dpca_norm(double *dtrj, double *sctrj, int numstep, int numdihed);
int dpca_covm(double *sctrj_n, int numstep, int numdihed, double *cov);
int dpca_diag(double *cov,double *eigenval,int numdihed);
void dpca_proj(double *sctraj_n,double *dpca,double *U,int numstep, int numdihed);
int dpca_proj_wdim(double *sctrj, double *dtrjv,double *vec ,int numstep, int numdihed,int numdim);
