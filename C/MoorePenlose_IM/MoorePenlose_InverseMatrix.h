
#ifndef INCLUDE_MP_inv
#define INCLUDE_MP_inv

int svd(double *mat, int m, int n, double *matU, double *matV, double *sv);

void mtrans(double *m1, double *m2, int n);

void mntrans(double *m1, double *m2, int m, int n);

void mnmult(double *mat1,int m1, int n1, double *mat2,int m2, int n2, double *m1m2);

void MPginvm(double *mat, double *invmat, int m, int n);

#endif
