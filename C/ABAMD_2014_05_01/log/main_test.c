#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "f2c.h"
#include "blaswrap.h"
#include "clapack.h"

int main(int argc, char *argv[])
{
  int i,j,k;
  double *mattemp,*test;
  static long int m,n,lda,info,piv[500],lwork=500;
  static double work[500];

  double *mat, *invmat;
  int num;

  //  DGETRF(&m,&n,mattemp,&lda,piv,&info);
  dgetrf_(&m,&n,mattemp,&lda,piv,&info);
  if (info!=0) return 0;
  //  DGETRF(&n,mattemp,&lda,piv,work,&lwork,&info);
  dgetri_(&n,mattemp,&lda,piv,work,&lwork,&info);
  if (info!=0) return 0;
   
  return 0;
}
