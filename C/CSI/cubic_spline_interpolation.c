
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "CSI.h"
#include "EF.h"
#include "IO.h"

#include "f2c.h"
#include "clapack.h"

//void mtrans(double *m1, double *m2, int n);
void csi_mtrans(double *m1,int n);

void csi_mmult(double *matrix,double *u,int N);

double csi_solv_matrix(double *matrix, double *v, double *u, int numpoint) {
  int i,j;

  //  double *matrixT;
  //  double *matrix2;
  long N,NRHS,lda,LDB,*IPIV,info;

  N     = numpoint-2;
  lda   = N;
  NRHS  = 1;
  LDB   = N;
  IPIV  = (long *)gcemalloc(sizeof(long)*N);
  // matrixT = (double *)gcemalloc(sizeof(double)*N*N);
  /*****************************************/
  /* for (i=0;i<N;++i)			   */
  /*   for (j=0;j<N;++j)		   */
  /*     matrix2[i*N+j]=matrix[i*N+j];	   */
  /* for (i=0;i<N;++i) {		   */
  /*   for (j=0;j<N;++j) {		   */
  /*     printf("%6.4lf ",matrix2[i*N+j]); */
  /*   }				   */
  /*   printf("\n");			   */
  /* }					   */
  /* 					   */
  /* for (i=0;i<N;++i) {		   */
  /*   printf("%6.4lf\n",v[i]);		   */
  /* }					   */
  /* printf("\n");			   */
  /*****************************************/
  csi_mtrans(matrix,N);
  dgesv_(&N,&NRHS,matrix/*T*/,&lda,IPIV,v,&LDB,&info);
  if (info!=0){
    printf("error; dgesv calculation!!\n info=%d\n",info);
    exit(1);
  }
  //  mmult(matrix2,v,N);

  u[0]=0.0;
  u[numpoint-1]=0.0;
  for (i=0;i<N;++i)
    u[i+1]=v[i];
}

double csi_set_matrix(double *x, double *y, double *matrix, double *v, int numpoint ) {
  int i,j,k;
  double *h;

  int nmat;

  h=(double *)gcemalloc(sizeof(double)*numpoint-1);

  for (i=0;i<numpoint-1;++i) h[i]=x[i+1]-x[i];
  for (i=1;i<numpoint-1;++i) v[i-1]=6.0*((y[i+1]-y[i])/h[i]-(y[i]-y[i-1])/h[i-1]);

  nmat=numpoint-2;
  for (i=0;i<nmat*nmat;++i) matrix[i]=0.0;
  for (i=0;i<nmat;++i) matrix[i*nmat+i]=2.0*(h[i]+h[i+1]);
  for (i=0;i<nmat;++i) {
    if ((i+1)<nmat) {
      matrix[i*nmat+(i+1)]=h[i+1];
      matrix[(i+1)*nmat+i]=h[i+1];
    }
  }
  
}

double csi_set_coff(double *u, double *x, double *y,  double *a, double *b, double *c, double *d, int numpoint ) {
  int i;

  for (i=0;i<numpoint-1;++i) {
    d[i]=y[i];
    b[i]=0.5*u[i];
  }

  for (i=0;i<numpoint-1;++i) {
    a[i]=(u[i+1]-u[i])/(6.0*(x[i+1]-x[i]));
    c[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(x[i+1]-x[i])/6.0*(2.0*u[i]+u[i+1]);
  }
}

double csi_calc_abcd(double *x, double *y, double *a, double *b, double *c, double *d, int numpoint ) {
  int i,j,k,l;
  double *matrix,*v,*u;

  matrix=(double*)gcemalloc(sizeof(double)*(numpoint-2)*(numpoint-2));
  v=(double*)gcemalloc(sizeof(double)*(numpoint-1));
  u=(double*)gcemalloc(sizeof(double)*numpoint);

  csi_set_matrix(x,y,matrix,v,numpoint);
  csi_solv_matrix(matrix,v,u,numpoint);
  csi_set_coff(u,x,y,a,b,c,d,numpoint);
}

double csi_interpolation(double *x, double *a, double *b, double *c, double *d,double *x_interpolated, double *y_interpolated,double width, int numpoint,int numinterpolated ) {
  int i,j;

  j=0;
  for (i=0;i<numinterpolated;++i) {
    x_interpolated[i]=width*i+x[0];
    
    for (;j<numinterpolated;++j)
      if ( (x_interpolated[i] >= x[j]) && (x_interpolated[i] < x[j+1]) )
	   break;

    y_interpolated[i]=a[j]*(x_interpolated[i]-x[j])*(x_interpolated[i]-x[j])*(x_interpolated[i]-x[j])
      +b[j]*(x_interpolated[i]-x[j])*(x_interpolated[i]-x[j])
      +c[j]*(x_interpolated[i]-x[j])
      +d[j];
  }

}

double csi_interpolation2(double *x, double *a, double *b, double *c, double *d,double *x_interpolated, double *y_interpolated,double width, int numpoint,int numinterpolated ) {
  int i,j;

  j=0;
  for (i=0;i<numinterpolated;++i) {
    for (;j<numinterpolated;++j)
      if ( (x_interpolated[i] >= x[j]) && (x_interpolated[i] <= x[j+1]) )
	   break;

    y_interpolated[i]=a[j]*(x_interpolated[i]-x[j])*(x_interpolated[i]-x[j])*(x_interpolated[i]-x[j])
      +b[j]*(x_interpolated[i]-x[j])*(x_interpolated[i]-x[j])
      +c[j]*(x_interpolated[i]-x[j])
      +d[j];
  }

}


void csi_mtrans(double *m1,int n){
  int i,j;
  double temp;

  for (i=0;i<n;++i) {
    for (j=i+1;j<n;++j) {
      temp = m1[i*n+j];
      m1[i*n+j] = m1[j*n+i];
      m1[j*n+i] = temp;
    }
  }
}

void csi_mmult(double *matrix,double *u,int N){
  int i,j;
  double *vec;

  vec=(double *)gcemalloc(sizeof(double)*N);

  for (j=0;j<N;++j)
    vec[i] = 0.0;

  for (i=0;i<N;++i)
    for (j=0;j<N;++j)
      vec[i] += matrix[i*N+j]*u[j];

  printf("\n");
  for (i=0;i<N;++i)
    printf("%6.4lf\n",vec[i]);

}

