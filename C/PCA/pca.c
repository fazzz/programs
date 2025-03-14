
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "EF.h"
#include "LA.h"
#include "PCA.h"
#include "f2c.h"
#include "clapack.h"

double kb=1.98723e-3*4.18407*100.0;

void mtrans(double *m1, double *m2, int n);

int pca_diag(double *cov,double *eigenval,int numatom) {
  int i,j;
  char jobz='V';
  char uplo='U';
  double *w,*work,*covm;
  long int n,lda,lwork,info;
  double var;

  n     = numatom*3;
  lda   = n;
  lwork = numatom*3*3;
  w     = (double *)gcemalloc(sizeof(double)*numatom*3);
  work  = (double *)gcemalloc(sizeof(double)*numatom*3*3);
  covm  = (double *)gcemalloc(sizeof(double)*numatom*3*numatom*3);
  mtrans(cov,covm,numatom*3);
  dsyev_(&jobz,&uplo,&n,covm,&lda,w,work,&lwork,&info);
  if (info!=0){
    printf("error; eigen vector calculation!!\n info=%d\n",info);
    exit(1);
  }
  for (i=0;i<numatom*3;++i)
    eigenval[numatom*3-1-i]=w[i];
  for (i=0;i<numatom*3;++i)
    for (j=0;j<numatom*3;++j)
      cov[i*numatom*3+j]=covm[i*numatom*3+j];

  return 1;
}

int pca_proj(double *traj,double *vec ,int numstep, int numatom) {
  int i,j,k;
  double *trajv;

  trajv=(double *)gcemalloc(sizeof(double)*numatom*3*numstep);
  for (i=0;i<numstep;++i)
    for (j=0;j<numatom*3;++j)
      for (k=0;k<numatom*3;++k)
	trajv[i*numatom*3+j]+=vec[j*numatom*3+k]*traj[i*numatom*3+k];

  for (i=0;i<numatom*3;++i)
    for (j=0;j<numatom*3;++j)
      traj[i*numatom*3+j]=trajv[i*numatom*3+j];

  return 1;
}

int pca_proj_wdim(double *traj,double *vec ,int numstep, int numatom,int numdim) {
  int i,j,k;
  double *trajv;

  trajv=(double *)gcemalloc(sizeof(double)*numatom*3*numstep);
  for (i=0;i<numstep;++i)
    for (j=0;j<numdim;++j)
      for (k=0;k<numatom*3;++k)
	trajv[i*numatom*3+j]+=vec[j*numatom*3+k]*traj[i*numatom*3+k];

  for (i=0;i<numatom*3;++i)
    for (j=0;j<numatom*3;++j)
      traj[i*numatom*3+j]=trajv[i*numatom*3+j];

  return 1;
}

int pca_covm(double *traj_a, int numstep, int numatom, double *cov){
  int i,j,k;

  for (i=0;i<numatom*3;++i)
    for (j=i;j<numatom*3;++j)
      cov[i*numatom*3+j]=0.0;
  for (i=0;i<numstep;++i)
    for (j=0;j<numatom*3;++j)
      for (k=j;k<numatom*3;++k)
	cov[j*numatom*3+k]=(i*cov[j*numatom*3+k]+traj_a[i*numatom*3+j]*traj_a[i*numatom*3+k])/(i+1);

  return 1;
}

int pca_norm(double *traj, int numstep, int numatom){
  int i,j;
  double *ave;

  ave  = (double *)ecalloc(sizeof(double),numatom*3);
  for (i=0;i<numstep;++i)
    for (j=0;j<numatom*3;++j)
      ave[j]=(i*ave[j]+traj[i*numatom*3+j])/(i+1);
  for (i=0;i<numstep;++i)
    for (j=0;j<numatom*3;++j)
      traj[i*numatom*3+j]=traj[i*numatom*3+j]-ave[j];

  return 0.0;
}

void pepca_norm(double *intene, int numstep, int numterm,double tp){
  int i,j;
  double *ave;

  ave=(double *)gcemalloc(sizeof(double)*numterm);
  for (i=0;i<numterm;++i)
    ave[i]=0.0;
  for (i=0;i<numstep;++i)
    for (j=0;j<numterm;++j)
      ave[j]=(i*ave[j]+intene[i*numterm+j])/(i+1);
  for (i=0;i<numstep;++i)
    for (j=0;j<numterm;++j)
      intene[i*numterm+j]=intene[i*numterm+j]-ave[j];

}

void pepca_covm(double *intene, int numstep, int numterm, double *cov){
  int i,j,k;

  for (i=0;i<numterm;++i)
    for (j=/*i*/0;j<numterm;++j)
      cov[i*numterm+j]=0.0;
  for (i=0;i<numstep;++i)
    for (j=0;j<numterm;++j)
      for (k=/*j*/0;k<numterm;++k)
	cov[j*numterm+k]=(i*cov[j*numterm+k]+intene[i*numterm+j]*intene[i*numterm+k])/(i+1);
}

void pepca_diag(double *cov,double *eigenval,int numterm){
  int i,j;
  char jobz='V';
  char uplo='U';
  double *w,*work,*covm;
  long int n,lda,lwork,info;
  double var;

  n     = numterm;
  lda   = n;
  lwork = numterm*3;
  w     = (double *)gcemalloc(sizeof(double)*numterm);
  work  = (double *)gcemalloc(sizeof(double)*numterm*3);
  covm  = (double *)gcemalloc(sizeof(double)*numterm*numterm);
  for (i=0;i<numterm;++i) {
    w[i]=0.0;work[i]=0.0;
    for (j=0;j<numterm;++j)
      covm[i*numterm+j]=0.0;
  }
  mtrans(cov,covm,numterm);
  dsyev_(&jobz,&uplo,&n,covm,&lda,w,work,&lwork,&info);
  if (info!=0){
    printf("error; eigen vector calculation!!\n info=%d\n",info);
    exit(1);
  }
  for (i=0;i<numterm;++i)
    eigenval[numterm-1-i]=w[i];
  for (i=0;i<numterm;++i)
    for (j=0;j<numterm;++j)
      cov[i*numterm+numterm-1-j]=covm[j*numterm+i];

}

void pepca_proj(double *intene,double *pepca,double *U,int numstep, int numterm) {
  int i,j,k;

  for (i=0;i<numstep;++i)
    for (j=0;j<numterm;++j)
      pepca[i*numterm+j]=0.0;
  for (i=0;i<numstep;++i)
    for (j=0;j<numterm;++j)
      for (k=0;k<numterm;++k)
	pepca[i*numterm+j]+=U[k*numterm+j]*intene[i*numterm+k];

}

void pepca_avevar(double *pepca,int numstep, int numterm,double *ave,double *var) {
  int i,j,k;

  for (i=0;i<numstep;++i) {
    for (j=0;j<numterm;++j) {
      ave[j]=(i*ave[j]+pepca[i*numterm+j])/(i+1);
      var[j]=(i*var[j]+pepca[i*numterm+j]*pepca[i*numterm+j])/(i+1);
    }
  }

  for (j=0;j<numterm;++j) {
    var[j]=var[j]-ave[j]*ave[j];
  }

}

void mtrans(double *m1, double *m2, int n){
  int i,j;

  for (i=0;i<n;++i)
    for (j=0;j<n;++j)
      m2[i*n+j] = m1[j*n+i];

}
