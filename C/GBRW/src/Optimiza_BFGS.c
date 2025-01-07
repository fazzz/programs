
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lbfgs.h"
#include "Optimiza_BFGS.h"
#include "EF.h"

/*static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    )
{
    int i;
    lbfgsfloatval_t fx = 0.0;
    for (i = 0;i < n;i += 2) {
        lbfgsfloatval_t t1 = 1.0 - x[i];
        lbfgsfloatval_t t2 = 10.0 * (x[i+1] - x[i] * x[i]);
        g[i+1] = 20.0 * t2;
        g[i] = -2.0 * (x[i] * g[i+1] + t1);
        fx += t1 * t1 + t2 * t2;
    }
    return fx;
    }*/

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *g_k/**w_k*/,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step)
{
  int i,j,k,l,ii,jj,kk;

  lbfgsfloatval_t f = 0.0;

  double pi,sqrt_pi;

  int K,N,*n_i;
  double *bk_i,h,*x_i,*x_k,**x_ij;
  double ***exp_kii_xij;

  double n1,n11,n2,n21,n3,n33,n4,n44,din,num;
  double *a_i,**b_ki,*c_k,*d_i,***e_ijk;

  struct data *dat = instance;

  K=(*dat).K;
  N=(*dat).n_sim;

  bk_i=(double *)gcemalloc(sizeof(double)*N);
  for (i=0;i<N;++i) bk_i[i]=(*dat).beta*(*dat).k[i];

  h=(*dat).h;

  n_i=(int *)gcemalloc(sizeof(int)*N);
  for (i=0;i<N;++i) n_i[i]=(*dat).n[i];

  x_ij=(double **)gcemalloc(sizeof(double *)*N);
  for (i=0;i<N;++i) x_ij[i]=(double *)gcemalloc(sizeof(double)*n_i[i]);

  for (i=0;i<N;++i) for (j=0;j<n_i[i];++j) x_ij[i][j]=(*dat).x_ij[i][j];

  x_i=(double *)gcemalloc(sizeof(double)*N);
  for (i=0;i<N;++i) x_i[i]=(*dat).x_i[i];

  x_k=(double *)gcemalloc(sizeof(double)*K);
  for (i=0;i<(*dat).K;++i) x_k[i]=(*dat).x_k[i];

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
  
  pi=acos(-1.0);

  a_i=(double *)gcemalloc(sizeof(double)*N);
  b_ki=(double **)gcemalloc(sizeof(double *)*K);
  for (k=0;k<K;++k) b_ki[k]=(double *)gcemalloc(sizeof(double)*N);
  c_k=(double *)gcemalloc(sizeof(double)*K);
  d_i=(double *)gcemalloc(sizeof(double)*N);
  e_ijk=(double ***)gcemalloc(sizeof(double **)*N);
  for (i=0;i<N;++i) {
    e_ijk[i]=(double **)gcemalloc(sizeof(double *)*n_i[i]);
    for (j=0;j<n_i[i];++j) e_ijk[i][j]=(double *)gcemalloc(sizeof(double)*K);
  }

  for (i=0;i<N;++i)  {
    a_i[i]=sqrt(2.0*pi/(1.0/h+bk_i[i]));
    d_i[i]=0.5*bk_i[i]*x_i[i]*x_i[i];
    for (k=0;k<K;++k) {
      b_ki[k][i]=0.5*(1.0/h*x_k[k]+bk_i[i]*x_i[i])*(1.0/h*x_k[k]+bk_i[i]*x_i[i])/(1.0/h+bk_i[i]);
    }
    for (j=0;j<n_i[i];++j) for (k=0;k<K;++k) e_ijk[i][j][k]=0.5*(x_ij[i][j]-x_k[k])*(x_ij[i][j]-x_k[k])/h;
  }
  for (k=0;k<K;++k)  c_k[k]=0.5/h*x_k[k]*x_k[k];
  
  n1=0.0;
  for (i=0;i<N;++i){
    n11=0.0; for (k=0;k<K;++k) n11+=/*sqrt(2.0*pi*h)**/a_i[i]*exp(g_k[k]-b_ki[k][i]-c_k[k]-d_i[i]);
    n1+=n_i[i]*log(n11);
  }

  n2=0.0;
  for (i=0;i<N;++i){
    for (j=0;j<n_i[i];++i){
      n21=0.0; for (k=0;k<K;++k) n21+=/*sqrt(2.0*pi*h)**/exp(g_k[k]-e_ijk[i][j][k]);
      n21=log(n21);
      n2+=n21;
    }
  }
  f=-n1+n2/*n1-n2*/;

  for (k=0;k<K;++k){
    n1=0.0;
    for (i=0;i<N;++i){
      num=exp(g_k[k]-b_ki[k][i]-c_k[k]);
      din=0.0; for (kk=0;kk<K;++kk) din+=exp(g_k[kk]-b_ki[kk][i]-c_k[kk]);
      n1+=n_i[i]*num/din;
    }

    n2=0.0;
    for (i=0;i<N;++i){
      for (j=0;j<n_i[i];++j){
	num=exp(g_k[k]-e_ijk[i][j][k]);
	din=0.0; for (kk=0;kk<K;++kk) din+=exp(g_k[kk]-e_ijk[i][j][kk]);
	n2+=num/din;
      }
    }
    g[k]=-n1+n2/*n1-n2*/;
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

int optimize_lnL_BFGS(double *g/**w*/, struct data dat ) {
  int i,j,k,ret = 0;
  lbfgsfloatval_t f;
  //  lbfgsfloatval_t *w_k = lbfgs_malloc(dat.K);
  lbfgsfloatval_t *g_k = lbfgs_malloc(dat.K);
  lbfgs_parameter_t param;

  if (g_k/*w_k*/ == NULL) {
    printf("ERROR: Failed to allocate a memory block for variables.\n");
    return 1;
  }

  for (i=0;i<dat.K;++i) g_k[i] = log(1.0/dat.K);
  
  lbfgs_parameter_init(&param);

  ret = lbfgs(dat.K, g_k/*w_k*/, &f, evaluate, progress,/*(void *) */&dat, &param);

  printf("L-BFGS optimization terminated with status code = %d\n", ret);

  //  for (i=0;i<dat.K;++i) w[i] = w_k[i];

  for (i=0;i<dat.K;++i) g[i] = g_k[i];

  //  lbfgs_free(w_k);

  lbfgs_free(g_k);

  return 0;
}
