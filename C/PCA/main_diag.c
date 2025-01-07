
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j,k;
  int flag=0;
  
  double *A,*diag,*temp,*test,*diagmat;

  A=(double *)egcmalloc(sizeof(double)*2*2);
  temp=(double *)egcmalloc(sizeof(double)*2*2); 
  test=(double *)egcmalloc(sizeof(double)*2*2);
  diag=(double *)egcmalloc(sizeof(double)*2);
  diagmat=(double *)egcmalloc(sizeof(double)*2*2);
  A[0*2+0]=111.34;
  A[0*2+1]=41.41;
  A[1*2+0]=41.41;
  A[1*2+1]=91.25;
  for (i=0;i<2;++i){
    for (j=0;j<2;++j){
      temp[i*2+j]=0.0;
      test[i*2+j]=0.0;
    }
  }
  for (i=0;i<2;++i)
    diag[i]=0.0;

  pepca_diag(A,diag,2);

  printf("eigenvalue=( ");
  for (i=0;i<2;++i)
    printf("%lf ",diag[i]);
  printf(") ");
  printf("\n\n");

  for (i=0;i<2;++i){
    for (j=0;j<2;++j){
      printf("%lf ",A[i*2+j]);
    }
    printf("\n");
  }

  for (i=0;i<2;++i){
    for (j=0;j<2;++j){
      if (i==j)
	diagmat[i*2+j]=diag[i];
      else
	diagmat[i*2+j]=0.0;
    }
  }
  printf("\n");

  /*****************/
  /* A[0]=-0.7860; */
  /* A[1]=0.6182;  */
  /* A[2]=-0.6182; */
  /* A[3]=-0.7860; */
  /*****************/

  diagmat[0]=111.34;
  diagmat[1]=41.41;
  diagmat[2]=41.41;
  diagmat[3]=91.25;

  for (i=0;i<2;++i) {
    for (j=0;j<2;++j) {
      temp[i*2+j]=0.0;
      test[i*2+j]=0.0;
    }
  }
  for (i=0;i<2;++i)
    for (j=0;j<2;++j)
      for (k=0;k<2;++k)
	temp[i*2+j]+=A[k*2+i]*diagmat[k*2+j];
  for (i=0;i<2;++i)
    for (j=0;j<2;++j)
      for (k=0;k<2;++k)
	test[i*2+j]+=temp[i*2+k]*A[k*2+j];

  for (i=0;i<2;++i){
    for (j=0;j<2;++j){
      printf("%lf ",test[i*2+j]);
    }
    printf("\n");
  }

  return 0;
}

