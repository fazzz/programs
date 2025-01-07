
void LAGC_v_product(double *v,double *mat);
void LAGC_mmult(double *m1, double *m2, double *m1m2, int n);
void LAGC_mnmult(double *mat1,int m1, int n1, double *mat2,int m2, int n2, double *m1m2);
void LAGC_mvmult(double *m, double *v, double *mv, int n);
void LAGC_mtrans(double *m1, double *m2, int n);
void LAGC_mntrans(double *m1, double *m2, int m, int n);
void LAGC_msetIni(double *m, int n);
void LAGC_msetzero(double *m, int n);
int LAGC_invm(double *mat, double *invmat, int num);
double LAGC_inprod(double *v1, double *v2, int n);
double LAGC_vtmvmult(double *vec1,double *mat,double *vec2,int num);
int LAGC_svd(double *mat, int m, int n, double *matU, double *matVT, double *sv);
int LAGC_MPginvm(double *mat, double *invmat, int m, int n);
int LAGC_MPginvm2(double *mat, double *invmat, int m, int n);
int LAGC_invm2(double *mat, double *invmat, int num);
