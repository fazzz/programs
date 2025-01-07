#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "f2c.h"
#include "clapack.h"

#include "MoorePenlose_InverseMatrix.h"

int svd(double *mat, int m, int n, double *matU, double *matV, double *sv) {
  int i,j;
  double *matT;
  double *U,*VT,*sigma,*Usigma,*A,*V,*UT;

  char jobu,jobvt;
  static long int m2,n2,lda,ldu,ldvt,info,lwork=1000/*200*/;
  static double work[1000/*200*/];
 
  jobu='a';
  jobvt='a';

  lda=m;
  ldu=m;

  ldvt=n;
  m2=m;
  n2=n;
  matT=(double *)malloc(sizeof(double)*m*n);
  mntrans(mat,matT,m,n);

  U=(double *)malloc(sizeof(double)*m*m);
  VT=(double *)malloc(sizeof(double)*n*n);

  dgesvd_(&jobu,&jobvt,&m2,&n2,matT,&lda,sv,U,&ldu,VT,&ldvt,work,&lwork,&info);

  if (info!=0) return 0;
  mtrans(U,matU,m);
  mtrans(VT,matV,n);
  V=malloc(sizeof(double)*n*n);
  for (i=0;i<n*n;++i) V[i]=VT[i];
  mtrans(V,VT,n);

  free(matT);
  free(U);
  free(VT);
  free(V);

}

void mtrans(double *m1, double *m2, int n){
  int i,j;

  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      m2[i*n+j] = m1[j*n+i];

}

void mntrans(double *m1, double *m2, int m, int n){
  int i,j,k;

  for (i=0;i<m*n;++i) {
    j=i%n;
    k=(int)(i/n);
    m2[j*m+k] = m1[i];
  }
}

void mnmult(double *mat1,int m1, int n1, double *mat2,int m2, int n2, double *m1m2){
  int i,j,k;

  if (n1!=m2) {
    printf("error\n");
    exit(1);
  }

  for (i=0;i<m1;++i)
    for (j=0;j<n2;++j)
      m1m2[i*n2+j] = 0.0;

  for (i=0;i<m1;++i)
    for (j=0;j<n2;++j)
      for (k=0;k<n1;++k)
	m1m2[i*n2+j] += mat1[i*n1+k]*mat2[k*n2+j];
}

void MPginvm(double *mat, double *invmat, int m, int n) {
  int i,j,k;
  int minmn;
  double *matU,*matVT,*matV,*matUT,*sv,*sigma,*VTsigma;
  double *AAP,*AAPA;

  matU=malloc(sizeof(double)*m*m);
  matUT=malloc(sizeof(double)*m*m);
  matVT=malloc(sizeof(double)*n*n);
  matV=malloc(sizeof(double)*n*n);

  minmn=n;  if (minmn > m) minmn=m;
  sv=malloc(sizeof(double)*minmn);

  svd(mat,m,n,matU,/*matV*/matVT,sv);

  ///////////////////////////////////////
  /**************************************************/
  /* printf("U= \n");				    */
  /* for (i=0;i<m;++i) {			    */
  /*   for (j=0;j<m;++j) {			    */
  /*     printf("%+5.3lf ",matU[i*m+j]);	    */
  /*   }      					    */
  /*   printf("\n");				    */
  /* }						    */
  /* printf("\n");				    */
  /* 						    */
  /* printf("VT= \n");				    */
  /* for (i=0;i<n;++i) {			    */
  /*   for (j=0;j<n;++j) {			    */
  /*     printf("%+5.3lf ",matVT/\*matV*\/[i*n+j]); */
  /*   }      					    */
  /*   printf("\n");				    */
  /* }						    */
  /* printf("\n");				    */
  /* 						    */
  /* printf("sv= \n");				    */
  /* for (i=0;i<minmn;++i) {			    */
  /*   printf("%+5.3lf ",sv[i]);		    */
  /* }						    */
  /* printf("\n");				    */
  /* printf("\n");				    */
  /**************************************************/
  ///////////////////////////////////////

  mtrans(matVT,matV,n);
  mtrans(matU,matUT,m);
  sigma=malloc(sizeof(double)*n*m);
  VTsigma=malloc(sizeof(double)*n*m);
  for (i=0;i<n*m;++i) sigma[i]=0.0;
  for (i=0;i<minmn;++i) {
    if (sv[i]>1.0e-10/*1.0e-7*/)
      sigma[i*n+i]=1.0/sv[i];
  }
  mnmult(matV,n,n,sigma,n,m,VTsigma);
  mnmult(VTsigma,n,m,matUT,m,m,invmat);

  free(matU);
  free(matUT);
  free(matVT);
  free(matV);
  free(sv);
  free(sigma);
  free(VTsigma);

}

/*
A =

1. 2. 3. 4. 
2. 4. 6. 8. 

V =
- 0.1825742 - 0.9689628 - 0.1609665 0.0432154 
- 0.3651484 0.1090576 - 0.4579506 - 0.8031528 
- 0.5477226 - 0.0405377 0.8070962 - 0.2166850 
- 0.7302967 0.2181152 - 0.3361053 0.5532863 

S =
12.247449 0. 0. 0. 
0. 9.730D-16 0. 0. 

U =
- 0.4472136 - 0.8944272 
- 0.8944272 0.4472136 

PinvA =
0.0066667 0.0133333 
0.0133333 0.0266667 
0.02 0.04 
0.0266667 0.0533333 

APA =
1. 2. 3. 4. 
2. 4. 6. 8. 

PAP =
0.0066667 0.0133333 
0.0133333 0.0266667 
0.02 0.04 
0.0266667 0.0533333
*/
