
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>

#include "MoorePenlose_InverseMatrix.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,m,n;

  double *A, *MPiA; 
  double *AMPiA,*AMPiAA;
  double *MPiAA,*MPiAAMPiA;

  m=2/*4*//*2*//*2*/;
  n=4/*4*//*3*//*2*/;

  A=(double *)malloc(sizeof(double)*m*n);
  MPiA=(double *)malloc(sizeof(double)*m*n);

  AMPiA=(double *)malloc(sizeof(double)*m*m);
  AMPiAA=(double *)malloc(sizeof(double)*m*n);
  MPiAA=(double *)malloc(sizeof(double)*n*n);
  MPiAAMPiA=(double *)malloc(sizeof(double)*n*m);

  A[0/*1,1*/]=1.0;  A[1/*1,2*/]=2.0;  A[2/*1,3*/]=3.0;  A[3/*1,4*/]=4.0;
  A[4/*2,1*/]=2.0;  A[5/*2,2*/]=4.0;  A[6/*2,3*/]=6.0;  A[7/*2,4*/]=8.0;

  //  A[0]=0.5;   A[1]=0.5;   A[2]=0.5;   A[3]=0.5;
  //  A[4]=0.5;   A[5]=0.5;   A[6]=-0.5;  A[7]=-0.5;
  //  A[8]=0.5;   A[9]=-0.5;  A[10]=0.5;  A[11]=-0.5;
  //  A[12]=0.5;  A[13]=-0.5; A[14]=-0.5; A[15]=0.5;

  /***************************/
  /* A[0]=1.0;   A[1]=1.0;   */
  /* A[2]=1.0;   A[3]=-1.0;  */
  /* A[4]=0.0;   A[5]=1.0;   */
  /***************************/

  /**************************/
  /* A[0]=1.0;   A[1]=1.0;  */
  /* A[2]=1.0;   A[3]=1.0;  */
  /**************************/

  printf("A= \n");
  for (i=0;i<m;++i) {
    for (j=0;j<n;++j) {
      printf("%+5.3lf ",A[i*n+j]);
    }      
    printf("\n");
  }
  printf("\n");

  MPginvm(A, MPiA, m, n);

  printf("A+= \n");
  for (i=0;i<n;++i) {
    for (j=0;j<m;++j) {
      printf("%+5.3lf ",MPiA[i*m+j]);
    }      
    printf("\n");
  }
  printf("\n");

  mnmult(A,m,n,MPiA,n,m,AMPiA);

  /****************************************/
  /* printf("AA+= \n");			  */
  /* for (i=0;i<m;++i) {		  */
  /*   for (j=0;j<m;++j) {		  */
  /*     printf("%+5.3lf ",AMPiA[i*m+j]); */
  /*   }      				  */
  /*   printf("\n");			  */
  /* }					  */
  /* printf("\n");			  */
  /****************************************/

  mnmult(AMPiA,m,m,A,m,n,AMPiAA);

  printf("AA+A= \n");
  for (i=0;i<m;++i) {
    for (j=0;j<n;++j) {
      printf("%+5.3lf ",AMPiAA[i*n+j]);
    }      
    printf("\n");
  }
  printf("\n");

  mnmult(MPiA,n,m,A,m,n,MPiAA);

  /****************************************/
  /* printf("A+A= \n");			  */
  /* for (i=0;i<n;++i) {		  */
  /*   for (j=0;j<n;++j) {		  */
  /*     printf("%+5.3lf ",MPiAA[i*n+j]); */
  /*   }      				  */
  /*   printf("\n");			  */
  /* }					  */
  /* printf("\n");			  */
  /****************************************/

  mnmult(MPiAA,n,n,MPiA,n,m,MPiAAMPiA);

  printf("A+AA+= \n");
  for (i=0;i<n;++i) {
    for (j=0;j<m;++j) {
      printf("%+5.3lf ",MPiAAMPiA[i*m+j]);
    }      
    printf("\n");
  }
  printf("\n");

  free(A);
  free(MPiA);
  
  return 0;
}
