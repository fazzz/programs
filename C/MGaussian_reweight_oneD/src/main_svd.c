#include <stdio.h>
#include "f2c.h"
#include "clapack.h"


#define ROW 6/*2*//*4*/

#define COL 5/*4*//*4*/

double A[ROW*COL]={      8.79,  6.11, -9.15,  9.57, -3.49,  9.84,
			 9.93,  6.91, -7.93,  1.64,  4.02,  0.15,
			 9.83,  5.04,  4.86,  8.83,  9.80, -8.99,
			 5.45, -0.27,  4.85,  0.74, 10.00, -6.02,
			 3.16,  7.98,  3.01,  5.80,  4.27, -5.31 };

int main(void) {

  static long int i,j;

  char jobu = 'A'/*'N'*/;

  char jobvt = 'A'/*'N'*/;

  static long int m = ROW, n = COL;

  static long int lda=ROW, ldu=/*1*/ROW, ldvt=/*1*/COL, lwork=100/*20*/, info;

  static double s[ROW], u[/*1*/ROW*ROW], vt[/*1*/COL*COL], work[100/*20*/];

  //  A[0]=0.5;A[1]=0.5;A[2]=0.5;A[3]=0.5;

  //  A[4]=0.5;A[5]=0.5;A[6]=-0.5;A[7]=-0.5;

  //  A[8]=0.5;A[9]=-0.5;A[10]=0.5;A[11]=-0.5;

  //  A[12]=0.5;A[13]=-0.5;A[14]=-0.5;A[15]=0.5;

  //  A[0]=1.0; A[1]=2.0; A[2]=3.0; A[3]=4.0;
  
  //  A[4]=2.0; A[5]=4.0; A[6]=6.0; A[7]=8.0;

  //  A[0]=1.0; A[2]=2.0; A[4]=3.0; A[6]=4.0;

  //  A[1]=2.0; A[3]=4.0; A[5]=6.0; A[7]=9.0;

  for(i=0;i<ROW;++i) {
    for(j=0;j<COL;++j) {
      printf("%+lf ",A[i*COL+j]);
    }
    printf("\n");
  }

  printf("\n");

  dgesvd_(&jobu, &jobvt, &m, &n, A, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);

  for(i=0;i<ROW;++i) printf("%lf\n",s[i]);

  printf("\n");

  for(i=0;i<COL;++i) {
    for(j=0;j<COL;++j) {
      printf("%+lf ",vt[i*COL+j]);
    }
    printf("\n");
  }

  printf("\n");

  for(i=0;i<ROW;++i) {
    for(j=0;j<ROW;++j) {
      printf("%+lf ",u[i*ROW+j]);
    }
    printf("\n");
  }

  return 0;

}

/*

matrix A:

     8.79   9.93   9.83   5.45   3.16
     6.11   6.91   5.04  -0.27   7.98
    -9.15  -7.93   4.86   4.85   3.01
     9.57   1.64   8.83   0.74   5.80
    -3.49   4.02   9.80  10.00   4.27
     9.84   0.15  -8.99  -6.02  -5.31


 Singular values
  27.47  22.64   8.56   5.99   2.01

 Left singular vectors (stored columnwise)
  -0.59   0.26   0.36   0.31   0.23
  -0.40   0.24  -0.22  -0.75  -0.36
  -0.03  -0.60  -0.45   0.23  -0.31
  -0.43   0.24  -0.69   0.33   0.16
  -0.47  -0.35   0.39   0.16  -0.52
   0.29   0.58  -0.02   0.38  -0.65

 Right singular vectors (stored rowwise)
  -0.25  -0.40  -0.69  -0.37  -0.41
   0.81   0.36  -0.25  -0.37  -0.10
  -0.26   0.70  -0.22   0.39  -0.49
   0.40  -0.45   0.25   0.43  -0.62
  -0.22   0.14   0.59  -0.63  -0.44
*/
