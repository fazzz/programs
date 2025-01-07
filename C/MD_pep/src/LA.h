#ifndef INCLUDE_LA
#define INCLUDE_LA

void v_product(double *v,double *mat);
double inprod(double *v1, double *v2, int n);
double outprod(double *v1,double *v2, double *v3);
void mvmult(double *m, double *v, double *mv, int n);

#endif
