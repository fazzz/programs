
void v_product(double *v,double *mat);
void LA_mmult(double *m1, double *m2, double *m1m2, int n);
void mnmult(double *mat1,int m1, int n1, double *mat2,int m2, int n2, double *m1m2);
void mvmult(double *m, double *v, double *mv, int n);
void mtrans(double *m1, double *m2, int n);
void mntrans(double *m1, double *m2, int m, int n);
void msetIni(double *m, int n);
void msetzero(double *m, int n);
int invm(double *mat, double *invmat, int num);
double inprod(double *v1, double *v2, int n);
double vtmvmult(double *vec1,double *mat,double *vec2,int num);
int svd(double *mat, int m, int n, double *matU, double *matVT, double *sv);
int MPginvm(double *mat, double *invmat, int m, int n);
int MPginvm2(double *mat, double *invmat, int m, int n);
int invm2(double *mat, double *invmat, int num);
int invm3(double *mat, double *invmat); // 2014-08-13
int invm4(double *a, double *inva); // 2014-08-13
void invm5(double *a, double *inv_a, int n); // 2014-08-13
