#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "LA.h"
#include "f2c.h"
#include "clapack.h"
#include "efunc.h"

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

  mat=(double *)ecalloc(sizeof(double),3*3);
  v_product(v1,mat);
  mvmult(mat,v2,v3,3);
  free(mat);
}

void mmult(double *m1, double *m2, double *m1m2, int n){
  int i,j,k;

  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      m1m2[i*n+j] = 0.0;

  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      for (k=0;k<n;++k)
	m1m2[i*n+j] += m1[i*n+k]*m2[k*n+j];
}

void mvmult(double *m, double *v, double *mv, int n){
  int i,j;

  for (i=0;i<n;++i)
    mv[i] = 0.0;
  
  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      mv[i] += m[i*n+j]*v[j];
}

void mtrans(double *m1, double *m2, int n){
  int i,j;

  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      m2[i*n+j] = m1[j*n+i];

}

void msetIni(double *m, int n) {
  int i;

  for (i=0;i<6;++i)
      m[i*n+i] = 1.0;
}

void msetzero(double *m, int n) {
  int i,j;

  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      m[i*n+j]=0.0;
}

int invm(double *mat, double *invmat, int num) {

  int i,j,k;
  double *mattemp,*test;
  static long int m,n,lda,info,piv[500],lwork=500;
  static double work[500];

  m = num;
  n = num;
  lda=num;
  
  mattemp=emalloc(sizeof(double)*m*n);
  test=emalloc(sizeof(double)*m*n);
  mtrans(mat,mattemp,m);
  dgetrf_(&m,&n,mattemp,&lda,piv,&info);
  if (info!=0) return 0;
  dgetri_(&n,mattemp,&lda,piv,work,&lwork,&info);
  if (info!=0) return 0;
  mtrans(mattemp,invmat,m);

  mmult(mat,invmat,test,3);
  
  for (i=0;i<3;++i){
    for (j=0;j<3;++j)
      printf("%lf ",test[i*3+j]);
    printf("\n");
  }

  free(mattemp);
  free(test);

  return 1;
}

double vtmvmult(double *vec1,double *mat,double *vec2,int num) {
  double *matvec2,vec1Tmatvec2;

  matvec2=(double *)emalloc(sizeof(double)*num);
  mvmult(mat,vec2,matvec2,num);
  vec1Tmatvec2=inprod(vec1,matvec2,num);

  return vec1Tmatvec2;
}
