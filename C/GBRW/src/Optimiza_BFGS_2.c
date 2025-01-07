
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lbfgs.h"
#include "Optimiza_BFGS_2.h"
#include "EF.h"

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *w_ij,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step)
{
  int i,j,k,ii,jj,iii,jjj;
  int l,ll,lll;

  lbfgsfloatval_t f = 0.0;

  double pi,sqrt_pi;

  int K,N,*n_i;
  double *bk_i,*x_i,**x_ij;
  double ***exp_kii_xij;

  double n1,n11,n2,n21,din,num;

  double **w_i_j,**g_ij;

  struct data *dat = instance;

  K=(*dat).K;
  N=(*dat).n_sim;

  bk_i=(double *)gcemalloc(sizeof(double)*N);
  for (i=0;i<N;++i) bk_i[i]=(*dat).beta*(*dat).k[i];

  n_i=(int *)gcemalloc(sizeof(int)*N);
  for (i=0;i<N;++i) n_i[i]=(*dat).n[i];

  x_ij=(double **)gcemalloc(sizeof(double *)*N);
  for (i=0;i<N;++i) x_ij[i]=(double *)gcemalloc(sizeof(double)*n_i[i]);

  for (i=0;i<N;++i) for (j=0;j<n_i[i];++j) x_ij[i][j]=(*dat).x_ij[i][j];

  x_i=(double *)gcemalloc(sizeof(double)*N);
  for (i=0;i<N;++i) x_i[i]=(*dat).x_i[i];

  exp_kii_xij=(double ***)gcemalloc(sizeof(double **)*N);
  for (i=0;i<N;++i) 
    exp_kii_xij[i]=(double **)gcemalloc(sizeof(double *)*N);
  for (i=0;i<N;++i) 
    for (j=0;j<N;++j) 
      exp_kii_xij[i][j]=(double *)gcemalloc(sizeof(double)*n_i[j]);

  for (i=0;i<N;++i)
    for (ii=0;ii<N;++ii)
      for (jj=0;jj<n_i[ii];++jj)
	exp_kii_xij[i][ii][jj]=(*dat).exp_kii_xij[i][ii][jj];

  w_i_j=(double **)gcemalloc(sizeof(double *)*N);
  for (i=0;i<N;++i) 
    w_i_j[i]=(double *)gcemalloc(sizeof(double)*n_i[i]);

  g_ij=(double **)gcemalloc(sizeof(double *)*N);
  for (i=0;i<N;++i) 
    g_ij[i]=(double *)gcemalloc(sizeof(double)*n_i[i]);

  l=0;
  for (i=0;i<N;++i) {
    for (j=0;j<n_i[i];++j) {
      w_i_j[i][j]=w_ij[l];
      ++l;
    }
  }
  
  pi=acos(-1.0);

  n1=0.0;
  n2=0.0;
  for (i=0;i<N;++i){
    n11=0.0;
    for (ii=0;ii<N;++ii){
      for (jj=0;jj<n_i[ii];++jj){
	n11+=exp(w_i_j[ii][jj])*exp_kii_xij[i][ii][jj];
      }
    }
    n1+=n_i[i]*log(n11);

    for (j=0;j<n_i[i];++j){
      //      n2+=log(w_i_j[i][j]);
      n2+=w_i_j[i][j];
    }
  }
  f=n1-n2;

  for (i=0;i<N;++i){
    for (j=0;j<n_i[i];++j){
      n1=0.0;
      for (ii=0;ii<N;++ii){
	num=exp_kii_xij[ii][i][j]*exp(w_i_j[i][j]);

	din=0.0;
	for (iii=0;iii<N;++iii) {
	  for (jjj=0;jjj<n_i[iii];++jjj) {
	    din+=exp(w_i_j[iii][jjj])*exp_kii_xij[ii][iii][jjj];
	  }
	}
	n1+=n_i[ii]*num/din;
      }
      //      n2=1.0/w_i_j[i][j];
      n2=1.0;
      g_ij[i][j]=n1-n2;
    }
  }

  l=0;
  for (i=0;i<N;++i) {
    for (j=0;j<n_i[i];++j) {
      g[l]=g_ij[i][j];
      ++l;
    }
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

  printf("Iteration %d:\n", k);
  printf("  fx = %f\n", fx);
  for (i=0;i<n;++i) printf("  x[%2d] = %f, \n",i+1, x[i]);
  printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
  printf("\n");

  return 0;
}

int optimize_lnL_BFGS(double *w, struct data dat ) {
  int i,j,k,ret = 0;
  lbfgsfloatval_t f;
  lbfgsfloatval_t *w_ij = lbfgs_malloc(dat.K);
  lbfgs_parameter_t param;

  if (w_ij == NULL) {
    printf("ERROR: Failed to allocate a memory block for variables.\n");
    return 1;
  }

  for (i=0;i<dat.K;++i) w_ij[i] = 1.0/dat.K;

  //  for (i=0;i<dat.K;++i) w_ij[i] = log(w_ij[i]);
  
  lbfgs_parameter_init(&param);

  ret = lbfgs(dat.K, w_ij, &f, evaluate, progress,/*(void *) */&dat, &param);

  printf("L-BFGS optimization terminated with status code = %d\n", ret);

  for (i=0;i<dat.K;++i) w[i] = w_ij[i];

  lbfgs_free(w_ij);

  return 0;
}
