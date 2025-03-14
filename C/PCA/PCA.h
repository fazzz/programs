
int pca_norm(double *traj, int numstep, int numatom);
int pca_covm(double *traj, int numstep, int numatom, double *cov);
int pca_diag(double *cov,double *eigenval,int numatom);
int pca_proj(double *traj,double *vec ,int numstep, int numatom);
int pca_proj_wdim(double *traj,double *vec ,int numstep, int numatom,int numdim);
void pepca_norm(double *intene, int numstep, int numterm,double tp);
void pepca_covm(double *intene, int numstep, int numatom, double *cov);
void pepca_diag(double *cov,double *eigenval,int numatom);
void pepca_proj(double *intene,double *pepca,double *U,int numstep, int numterm);
void pepca_avevar(double *pepca,int numstep, int numterm,double *ave,double *var);
