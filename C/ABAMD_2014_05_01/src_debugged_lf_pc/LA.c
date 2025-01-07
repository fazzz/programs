#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "LA.h"
#include "f2c.h"
#include "clapack.h"
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
  //  mat=(double *)gcemalloc(sizeof(double)*3*3); // 2014-07-22
  //  mat=(double *)emalloc(sizeof(double)*3*3); // 2014-07-22 // 2014-09-05
  mat=(double *)calloc(3*3,sizeof(double)); // 2014-07-22 // 2014-09-05
  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      mat[i*3+j]=0.0;
  v_product(v1,mat);
  mvmult(mat,v2,v3,3);
  //  free(mat);

  free(mat); // 2014-07-22

  return 0.0; // 2014-07-04
}

void LA_mmult(double *m1, double *m2, double *m1m2, int n){
  int i,j,k;

  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      m1m2[i*n+j] = 0.0;

  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      for (k=0;k<n;++k)
	m1m2[i*n+j] += m1[i*n+k]*m2[k*n+j];
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

void mntrans(double *m1, double *m2, int m, int n){
  int i,j,k;

  for (i=0;i<m*n;++i) {
    j=i%n;
    k=(int)(i/n);
    m2[j*m+k] = m1[i];
  }
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
  
  //  mattemp=emalloc(sizeof(double)*m*n);
  //  test=emalloc(sizeof(double)*m*n);
  //  mattemp=gcemalloc(sizeof(double)*m*n); // 2014-08-13
  //  test=gcemalloc(sizeof(double)*m*n);    // 2014-08-13
  mattemp=emalloc(sizeof(double)*m*n); // 2014-08-13

  mtrans(mat,mattemp,m);
  dgetrf_(&m,&n,mattemp,&lda,piv,&info);
  if (info!=0) return 0;
  dgetri_(&n,mattemp,&lda,piv,work,&lwork,&info);
  if (info!=0) return 0;
  mtrans(mattemp,invmat,m);

  //  mmult(mat,invmat,test,3);

  free(mattemp); // 2014-08-13
  //  free(mattemp);
  //  free(test);

  return 1;
}

int invm2(double *mat, double *invmat, int num) {

  int i,j,k;
  double *mattemp,*test;
  /******************************************************************/
  /* static long int m,n,lda,info,piv[500],lwork=500; // 2014-08-13 */
  /* static double work[500];                         // 2014-08-13 */
  /******************************************************************/
  static long int m,n,lda,info,piv[10],lwork=10;  // 2014-08-13
  static double work[10];                         // 2014-08-13
  /**********************************************************************/
  /* static long int m,n,lda,info,*piv,lwork=500;         // 2014-08-13 */
  /* static double *work;                                 // 2014-08-13 */
  /**********************************************************************/

  m = num;
  n = num;
  lda=num;

  /**************************************************************************/
  /* piv=(long int *)emalloc(sizeof(long int)*500); // 2014-08-13	    */
  /* work=(double *)emalloc(sizeof(double)*500); // 2014-08-13		    */
  /* info=0;  for (i=0;i<500;++i) { piv[i]=00; work[i]=0.0; } // 2014-08-13 */
  /**************************************************************************/
  
  //  mattemp=(double *)gcemalloc(sizeof(double)*m*n); // 2014-07-22
  //  mattemp=(double *)emalloc(sizeof(double)*m*n); // 2014-07-22 // 2014-09-05
  mattemp=(double *)calloc(m*n,sizeof(double)); // 2014-07-22 // 2014-09-05
  //  test=emalloc(sizeof(double)*m*n);
  mtrans(mat,mattemp,m);
  dgetrf_(&m,&n,mattemp,&lda,piv,&info);
  if (info!=0) return 0;
  dgetri_(&n,mattemp,&lda,piv,work,&lwork,&info);
  if (info!=0) return 0;
  mtrans(mattemp,invmat,m);

  //  mmult(mat,invmat,test,3);

  //  free(mattemp);
  //  free(test);

  free(mattemp); // 2014-07-22

  return 1;
}


double vtmvmult(double *vec1,double *mat,double *vec2,int num) {
  double *matvec2,vec1Tmatvec2;

  //  matvec2=(double *)emalloc(sizeof(double)*num); // 11_04_11
  //  matvec2=(double *)gcemalloc(sizeof(double)*num);   // 11_04_11 // 2014-07-22
  //  matvec2=(double *)emalloc(sizeof(double)*num);   // 2014-07-22 // 2014-09-05
  matvec2=(double *)calloc(num,sizeof(double));   // 2014-07-22 // 2014-09-05
  mvmult(mat,vec2,matvec2,num);
  vec1Tmatvec2=inprod(vec1,matvec2,num);

  free(matvec2); // 2014-07-22

  return vec1Tmatvec2;
}

int svd(double *mat, int m, int n, double *matU, double *matV, double *sv) {
  int i,j;
  double *matT;
  double *U,*VT,*sigma,*Usigma,*A,*V,*UT;
  //  double *A1,*A2,*A3,*A4;
  //  double sum=0.0;

  char jobu,jobvt;
  static long int m2,n2,lda,ldu,ldvt,info,lwork=10000;
  static double work[10000];
 
  jobu='A';
  jobvt='A';
  lda=m;
  ldu=m;
  ldvt=n;
  m2=m;
  n2=n;
  matT=(double *)gcemalloc(sizeof(double)*m*n);
  mntrans(mat,matT,m,n);

  U=(double *)gcemalloc(sizeof(double)*m*m);
  VT=(double *)gcemalloc(sizeof(double)*n*n);
  dgesvd_(&jobu,&jobvt,&m2,&n2,matT,&lda,sv,U,&ldu,VT,&ldvt,work,&lwork,&info);
  if (info!=0) return 0;
  mtrans(U,matU,m);
  mtrans(VT,matV,n);
  V=gcemalloc(sizeof(double)*n*n);
  for (i=0;i<n*n;++i) V[i]=VT[i];
  mtrans(V,VT,n);

  /**********************************************************************/
  /* sigma=gcemalloc(sizeof(double)*n*m);			        */
  /* Usigma=gcemalloc(sizeof(double)*n*m);			        */
  /* UT=gcemalloc(sizeof(double)*m*m);				        */
  /* mtrans(U,UT,m);						        */
  /* A=gcemalloc(sizeof(double)*m*n);				        */
  /* A1=gcemalloc(sizeof(double)*m*n);				        */
  /* A2=gcemalloc(sizeof(double)*m*n);				        */
  /* A3=gcemalloc(sizeof(double)*m*n);				        */
  /* A4=gcemalloc(sizeof(double)*m*n);				        */
  /* for (i=0;i<n*m;++i) sigma[i]=0.0;				        */
  /* for (i=0;i<n;++i) {					        */
  /*   if (sv[i]!=0.0)						        */
  /*     sigma[i*n+i]=sv[i];					        */
  /* }								        */
  /* mnmult(matU,m,m,sigma,m,n,Usigma);				        */
  /* mnmult(Usigma,m,n,VT,n,n,A1);				        */
  /* sum=0.0;for (i=0;i<m*n;++i) sum+=A1[i]-mat[i]; printf("%e\n",sum); */
  /**********************************************************************/
  /**********************************************************************/
  /* mnmult(UT,m,m,sigma,m,n,Usigma);				        */
  /* mnmult(Usigma,m,n,VT,n,n,A2);				        */
  /* sum=0.0;for (i=0;i<m*n;++i) sum+=A2[i]-mat[i]; printf("%e\n",sum); */
  /* mnmult(U,m,m,sigma,m,n,Usigma);				        */
  /* mnmult(Usigma,m,n,V,n,n,A3);				        */
  /* sum=0.0;for (i=0;i<m*n;++i) sum+=A3[i]-mat[i]; printf("%e\n",sum); */
  /* mnmult(UT,m,m,sigma,m,n,Usigma);				        */
  /* mnmult(Usigma,m,n,V,n,n,A4);				        */
  /* sum=0.0;for (i=0;i<m*n;++i) sum+=A4[i]-mat[i]; printf("%e\n",sum); */
  /**********************************************************************/
  /*****************************************/
  /************************************/
  /* mnmult(VT,n,n,sigma,n,m,Usigma); */
  /* mnmult(Usigma,n,m,VT,m,m,A);     */
  /************************************/

  /************************************/
  /* printf("USIVT=\n");  	      */
  /* for (i=0;i<m;++i) {	      */
  /*   for (j=0;j<n;++j) {	      */
  /*     printf("%10.4lf ",A[i*n+j]); */
  /*   }			      */
  /*   printf("\n");		      */
  /* }				      */
  /************************************/

  return 0; // 2014-07-04
}

int MPginvm(double *mat, double *invmat, int m, int n) {
  int i,j,k;
  int minmn;
  double *matU,*matVT,*matV,*matUT,*sv,*sigma,*VTsigma;
  double *AAP,*AAPA;

  matU=gcemalloc(sizeof(double)*m*m);
  matUT=gcemalloc(sizeof(double)*m*m);
  matVT=gcemalloc(sizeof(double)*n*n);
  matV=gcemalloc(sizeof(double)*n*n);
  minmn=n;
  if (minmn > m) minmn=m;
  sv=gcemalloc(sizeof(double)*minmn);

  svd(mat,m,n,matU,matV,sv);

  mtrans(matV,matVT,n);
  mtrans(matU,matUT,m);
  sigma=gcemalloc(sizeof(double)*n*m);
  VTsigma=gcemalloc(sizeof(double)*n*m);
  for (i=0;i<n*m;++i) sigma[i]=0.0;
  for (i=0;i<minmn;++i) {
    if (sv[i]>1.0e-4)
      sigma[i*n+i]=1.0/sv[i];
  }
  mnmult(matVT,n,n,sigma,n,m,VTsigma);
  mnmult(VTsigma,n,m,matU,m,m,invmat);

  // check for GIM calc.
  AAP=gcemalloc(sizeof(double)*m*m);
  AAPA=gcemalloc(sizeof(double)*m*n);
  
  mnmult(mat,m,n,invmat,n,m,AAP);
  mnmult(AAP,m,m,mat,m,n,AAPA);

  /**************************************/
  /* for (i=0;i<m;++i) {	        */
  /*   printf(":: ");		        */
  /*   for (j=0;j<n;++j)	        */
  /*     printf("%3.2lf ",mat[i*m+j]);  */
  /*   printf("\n ");		        */
  /*   printf("** ");		        */
  /*   for (j=0;j<n;++j)	        */
  /*     printf("%3.2lf ",AAPA[i*m+j]); */
  /*   printf("\n ");		        */
  /* }				        */
  /**************************************/

  return 0; // 2014-07-04
}

int MPginvm2(double *mat, double *invmat, int m, int n) {
  int i,j,k;
  int minmn;
  double *matU,*matVT,*matV,*matUT,*sv,*sigma,*VTsigma;
  double *AAP,*AAPA;

  //  double sum=0.0;

  matU=gcemalloc(sizeof(double)*m*m);
  matUT=gcemalloc(sizeof(double)*m*m);
  matVT=gcemalloc(sizeof(double)*n*n);
  matV=gcemalloc(sizeof(double)*n*n);
  minmn=n;
  if (minmn > m) minmn=m;
  sv=gcemalloc(sizeof(double)*minmn);

  svd(mat,m,n,matU,matVT,sv);

  mtrans(matVT,matV,n);
  mtrans(matU,matUT,m);
  sigma=gcemalloc(sizeof(double)*n*m);
  VTsigma=gcemalloc(sizeof(double)*n*m);
  for (i=0;i<n*m;++i) sigma[i]=0.0;
  for (i=0;i<minmn;++i) {
    if (sv[i]>1.0e-7)
      sigma[i*n+i]=1.0/sv[i];
  }
  mnmult(matV,n,n,sigma,n,m,VTsigma);
  mnmult(VTsigma,n,m,matUT,m,m,invmat);

  // check for GIM calc.
  /***************************************/
  /* AAP=gcemalloc(sizeof(double)*m*m);	 */
  /* AAPA=gcemalloc(sizeof(double)*m*n); */
  /* 					 */
  /* mnmult(mat,m,n,invmat,n,m,AAP);	 */
  /* mnmult(AAP,m,m,mat,m,n,AAPA);	 */
  /***************************************/
  
  //  sum=0.0;for (i=0;i<m*n;++i) sum+=AAPA[i]-mat[i]; printf("%e\n",sum);

  /*************************************/
  /* for (i=0;i<m;++i) {	       */
  /*   printf(":: ");		       */
  /*   for (j=0;j<n;++j)	       */
  /*     printf("%3.2e ",mat[i*m+j]);  */
  /*   printf("\n ");		       */
  /*   printf("** ");		       */
  /*   for (j=0;j<n;++j)	       */
  /*     printf("%3.2e ",AAPA[i*m+j]); */
  /*   printf("\n ");		       */
  /* }				       */
  /*************************************/

  return 0; // 2014-07-04
}

// 2014-08-13
int invm3(double *mat, double *invmat) {

  int i,j,k;
  double mattemp[9];
  static long int m,n,lda,info,piv[3],lwork=3;
  static double work[3];

  m = 3;
  n = 3;
  lda=3;

  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      mattemp[i*3+j] = mat[j*3+i];
  dgetrf_(&m,&n,mattemp,&lda,piv,&info);
  if (info!=0) return 0;
  dgetri_(&n,mattemp,&lda,piv,work,&lwork,&info);
  if (info!=0) return 0;
  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      invmat[i*3+j] = mattemp[j*3+i];

  return 1;
}

// 2014-08-13
int invm4(double *a, double *inva) {

  double det;

  det=a[0*3+0]*a[1*3+1]*a[2*3+2];
  det+=a[1*3+0]*a[2*3+1]*a[0*3+2];
  det+=a[2*3+0]*a[0*3+1]*a[1*3+2];
  det-=a[2*3+0]*a[1*3+1]*a[0*3+2];
  det-=a[1*3+0]*a[0*3+1]*a[2*3+2];
  det-=a[0*3+0]*a[2*3+1]*a[1*3+2];

  inva[0*3+0] = (a[1*3+1]*a[2*3+2]-a[1*3+2]*a[2*3+1])/det;
  inva[0*3+1] = (a[0*3+2]*a[2*3+1]-a[0*3+1]*a[2*3+2])/det;
  inva[0*3+2] = (a[0*3+1]*a[1*3+2]-a[0*3+2]*a[1*3+1])/det;

  inva[1*3+0] = (a[1*3+2]*a[2*3+0]-a[1*3+0]*a[2*3+2])/det;
  inva[1*3+1] = (a[0*3+0]*a[2*3+2]-a[0*3+2]*a[2*3+0])/det;
  inva[1*3+2] = (a[0*3+2]*a[1*3+0]-a[0*3+0]*a[1*3+2])/det;

  inva[2*3+0] = (a[1*3+0]*a[2*3+1]-a[1*3+1]*a[2*3+0])/det;
  inva[2*3+1] = (a[0*3+1]*a[2*3+0]-a[0*3+0]*a[2*3+1])/det;
  inva[2*3+2] = (a[0*3+0]*a[1*3+1]-a[0*3+1]*a[1*3+0])/det;

  return 1;
}

void invm5(double *a, double *inv_a, int n) {
  //  double a[4][4]={{1,2,0,-1},{-1,1,2,0},{2,0,1,1},{1,-2,-1,1}};
  //  double inv_a[4][4]; 
  double buf; 
  int i,j,k; 
  //  int n=4;  

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      inv_a[i*n+j]=(i==j)?1.0:0.0;
    }
  }

  for(i=0;i<n;i++){
    buf=1/a[i*n+i];
    for(j=0;j<n;j++){
      a[i*n+j]*=buf;
      inv_a[i*n+j]*=buf;
    }
    for(j=0;j<n;j++){
      if(i!=j){
	buf=a[j*n+i];
	for(k=0;k<n;k++){
	  a[j*n+k]-=a[i*n+k]*buf;
	  inv_a[j*n+k]-=inv_a[i*n+k]*buf;
	}
      }
    }
  }

/**********************************/
/* for(i=0;i<n;i++){		  */
/*   for(j=0;j<n;j++){		  */
/*     printf(" %f",inv_a[i][j]); */
/*   }				  */
/*   printf("\n");		  */
/**********************************/
}

/* 
 2.000000 2.000000 -1.000000 3.000000
 -4.000000 -5.000000 3.000000 -7.000000
 3.000000 4.000000 -2.000000 5.000000
 -7.000000 -8.000000 5.000000 -11.000000
*/
