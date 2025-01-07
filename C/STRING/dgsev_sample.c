
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "IO.h"

void mtrans(double *m1, double *m2, int n);
double csi_solv_matrix(double *matrix, double *v, double *u, int numpoint );

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int flag;

  int numpoint=3;
  double *matrix,*v,*u;

  matrix = (double *)gcemalloc(sizeof(double)*numpoint*numpoint);
  v = (double *)gcemalloc(sizeof(double)*numpoint);
  u = (double *)gcemalloc(sizeof(double)*numpoint);

  matrix[0]=1.0;  matrix[1]=1.0;  matrix[2]=1.0;
  matrix[3]=3.0;  matrix[4]=1.0;  matrix[5]=-3.0;
  matrix[6]=1.0;  matrix[7]=-2.0;  matrix[8]=-5.0;

  v[0]=3.0;
  v[1]=1.0;
  v[2]=-6.0;

  csi_solv_matrix(matrix,v,u,numpoint);

  for (i=0;i<numpoint;++i) 
    printf("%10.4lf n",u[i]);
}

double csi_solv_matrix(double *matrix, double *v, double *u, int numpoint ) {
  int i;

  double *matrixT;
  long N,NRHS,lda,LDB,*IPIV,info;

  N     = numpoint;
  lda   = N;
  NRHS  = N;
  LDB   = N;
  IPIV  = (long *)gcemalloc(sizeof(long)*N);
  matrixT = (double *)gcemalloc(sizeof(double)*N*N);
  mtrans(matrix,matrixT,N);
  dgesv_(&N,&NRHS,matrixT,&lda,IPIV,v,&LDB,&info);
  if (info!=0){
    printf("error; dgesv calculation!!\n info=%d\n",info);
    exit(1);
  }

  for (i=0;i<N;++i)
    u[i]=v[i];
}

void mtrans(double *m1, double *m2, int n){
  int i,j;

  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      m2[i*n+j] = m1[j*n+i];

}
