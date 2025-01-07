#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "LA.h"
#include "EF.h"

void v_product(double *v,double *mat) {
  mat[0]=  0.0;
  mat[1]=-v[2];
  mat[2]= v[1];
  mat[3]= v[2];
  mat[4]=  0.0;
  mat[5]=-v[0];
  mat[6]=-v[1];
  mat[7]= v[0];
  mat[8]=  0.0;
}

double inprod(double *v1, double *v2, int n) {
  int i;
  double in=0.0;

  for (i=0;i<n;++i)
    in += v1[i]*v2[i];

  return in;
}

double outprod(double *v1,double *v2, double *v3) {
  double *mat;
  int i,j;

  //  mat=(double *)ecalloc(sizeof(double),3*3);
  mat=(double *)gcemalloc(sizeof(double)*3*3);
  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      mat[i*3+j]=0.0;
  v_product(v1,mat);
  mvmult(mat,v2,v3,3);
  //  free(mat);
}

void mvmult(double *m, double *v, double *mv, int n){
  int i,j;

  for (i=0;i<n;++i)
    mv[i] = 0.0;
  
  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      mv[i] += m[i*n+j]*v[j];
}

