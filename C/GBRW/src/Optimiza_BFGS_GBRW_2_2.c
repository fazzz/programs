
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lbfgs.h"
#include "Optimiza_BFGS_GBRW_2.h"
#include "EF.h"

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *g_k,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step)
{
  int i,j,k,l,ii,jj,kk;

  lbfgsfloatval_t f = 0.0;

  double n1,n2,n11,n22,num,din;

  double pi,sqrt_pi;

  int K,N,*n_i;
  double h,*x_i,**x_ij;
  double *A_i,**B_ik,***S_k_ij;

  struct data *dat = instance;

  K=(*dat).K;
  N=(*dat).n_sim;
  h=(*dat).h;

  n_i=(int *)gcemalloc(sizeof(int)*N);
  for (i=0;i<N;++i) n_i[i]=(*dat).n[i];

  x_ij=(double **)gcemalloc(sizeof(double *)*N);
  for (i=0;i<N;++i) x_ij[i]=(double *)gcemalloc(sizeof(double)*n_i[i]);
  for (i=0;i<N;++i) for (j=0;j<n_i[i];++j) x_ij[i][j]=(*dat).x_ij[i][j];

  A_i=(double *)gcemalloc(sizeof(double)*N);
  for (i=0;i<N;++i) A_i[i]=(*dat).A_i[i];

  B_ik=(double **)gcemalloc(sizeof(double *)*N);
  for (i=0;i<N;++i) B_ik[i]=(double *)gcemalloc(sizeof(double *)*K);
  for (i=0;i<N;++i) for (j=0;j<K;++j) B_ik[i][j]=(*dat).B_ik[i][j];

  S_k_ij=(double ***)gcemalloc(sizeof(double **)*K);
  for (i=0;i<K;++i) S_k_ij[i]=(double **)gcemalloc(sizeof(double *)*N);
  for (i=0;i<K;++i) for (j=0;j<N;++j) S_k_ij[i][j]=(double *)gcemalloc(sizeof(double)*n_i[j]);
  for (k=0;k<K;++k)
    for (i=0;i<N;++i)
      for (j=0;j<n_i[i];++j)
	S_k_ij[k][i][j]=(*dat).S_k_ij[k][i][j];
    
  pi=acos(-1.0);

  f=0.0;
  n1=0.0;
  n2=0.0;
  for (i=0;i<N;++i){
    n11=0.0; for (k=0;k<K;++k) n11+=exp(g_k[k])*A[k]*C_ik[i][k]; n11=log(n11);
    n1+=n_i[i]*n11;

    for (j=0;j<n_i[i];++j){
      n22=0.0; for (k=0;k<K;++k) n22+=exp(g_k[k])*A[k]*S_k_ij[k][i][j]; n22=log(n22);
      n2+=n22;
    }
  }
  f=n1-n2;

  for (k=0;k<K;++k){
    n1=0.0;
    for (i=0;i<N;++i){
      num=exp(g_k[k])*A[k]*C_ik[i][k];
      din=0.0; for (kk=0;kk<K;++kk) din+=exp(g_k[kk])*A[kk]*C_ik[i][kk];
      n1+=n_i[i]*num/din;
    }

    n2=0.0;
    for (i=0;i<N;++i){
      for (j=0;j<n_i[i];++j){
	num=exp(g_k[k])*A[k]*S_k_ij[k][i][j];
	din=0.0; for (kk=0;kk<K;++kk) din+=exp(g_k[kk])*A[kk]*S_k_ij[kk][i][j];
	n2+=num/din;
      }
    }
    g[k]=n1-n2;
  }

  return f;
}

static int progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    )
{
  int i;

  /************************************************************************/
  /* printf("Iteration %d:\n", k);					  */
  /* printf("  fx = %f\n", fx);						  */
  /* for (i=0;i<n;++i) printf("  x[%2d] = %f, \n",i+1, x[i]);		  */
  /* printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step); */
  /* printf("\n");							  */
  /************************************************************************/

  printf("Iteration %d:  ", k);
  printf("  fx = %f\n", fx);

  return 0;
}

int optimize_lnL_BFGS_2(double *g, struct data dat ) {
  int i,j,k,ret = 0;
  lbfgsfloatval_t f;
  lbfgsfloatval_t *g_k = lbfgs_malloc(dat.K);
  lbfgs_parameter_t param;

  if (g_k == NULL) {
    printf("ERROR: Failed to allocate a memory block for variables.\n");
    return 1;
  }

  for (i=0;i<dat.K;++i) g_k[i] = log(1.0/dat.K);
  
  lbfgs_parameter_init(&param);

  ret = lbfgs(dat.K, g_k, &f, evaluate, progress,&dat, &param);

  printf("L-BFGS optimization terminated with status code = %d\n", ret);

  for (i=0;i<dat.K;++i) g[i] = g_k[i];

  lbfgs_free(g_k);

  return 0;
}
