
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "f2c.h"
#include "clapack.h"

int main(int argc, char *argv[]) {
  int i,j;
  double *Inertia,*InvInertia;

  Inertia=(double *)calloc(sizeof(double),9);
  InvInertia=(double *)calloc(sizeof(double),9);

  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      Inertia[i*3+j]=0.0;
  Inertia[0]=1.0;
  Inertia[4]=1.0;
  Inertia[8]=1.0;
  invm2(Inertia,InvInertia,3);


}

int invm2(double *mat, double *invmat, const int num) {

  int i,j,k;
  double *mattemp;
  static long int m1,n1,lda,info,piv[500],lwork=500;
  static double work[500];

  m1 = num;
  n1 = num;
  lda=m1;
  
  mattemp=malloc(sizeof(double)*m1*n1);
  memcpy(mattemp,mat,sizeof(double)*m1*n1);
  dgetrf_(&m1,&n1,mattemp,&lda,piv,&info);
  if (info!=0) return 0;
  dgetri_(&n1,mattemp,&lda,piv,work,&lwork,&info);
  if (info!=0) return 0;
  memcpy(mat,mattemp,sizeof(double)*m1*n1);

  return 1;
}
