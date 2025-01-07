
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lbfgs.h"
#include "Optimiza_BFGS_2.h"
#include "EF.h"

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *g_i,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step)
{
  int i,j,ii,jj,iii;
  int k,l;

  lbfgsfloatval_t f = 0.0;

  double pi,sqrt_pi;

  int K,N,*n_i;
  double *bk_i,*x_i,**x_ij;
  double ***exp_kii_xij;

  double n1,n11,n2,n21,din,num;

  double *f_i;
  double **P_ij,sum,lnp;

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
  
  pi=acos(-1.0);

  n1=0.0; for (i=0;i<N;++i) n1+=n_i[i]*g_i[i];
  
  n2=0.0;
  for (i=0;i<N;++i) {
    for (j=0;j<n_i[i];++j){
      n11=0.0; for (ii=0;ii<N;++ii) n11+=n_i[ii]*exp(g_i[ii])*exp_kii_xij[ii][i][j];
      n2+=log(n11);
    }
  }
  
  f=-n1+n2;

  /**********************************************************************************/
  /* f_i=(double *)gcemalloc(sizeof(double)*N);					    */
  /* for (i=0;i<N;++i) f_i[i]=exp(g_i[i]);					    */
  /* 										    */
  /* P_ij=(double **)gcemalloc(sizeof(double)*N);				    */
  /* for (i=0;i<N;++i) P_ij[i]=(double *)gcemalloc(sizeof(double)*n_i[i]);	    */
  /* 										    */
  /* for (i=0;i<N;++i) {							    */
  /*   for (j=0;j<n_i[i];++j) {							    */
  /*     P_ij[i][j]=0.0;							    */
  /*     for (k=0;k<N;++k) {							    */
  /* 	P_ij[i][j]+=n_i[k]*f_i[k]*exp_kii_xij[k][i][j];				    */
  /*     }									    */
  /*     P_ij[i][j]=1.0/P_ij[i][j];						    */
  /*   }									    */
  /* }										    */
  /* 										    */
  /* //  sum=0.0; for (i=0;i<N;++i) for (j=0;j<n_i[i];++j) sum+=P_ij[i][j];	    */
  /* //  for (i=0;i<N;++i) for (j=0;j<n_i[i];++j) P_ij[i][j]=P_ij[i][j]/sum;	    */
  /* 										    */
  /* f=0.0;									    */
  /* for (i=0;i<N;++i) {							    */
  /*   for (j=0;j<n_i[i];++j) {							    */
  /*     lnp=f_i[i]*exp_kii_xij[i][i][j]*P_ij[i][j];				    */
  /*     lnp=log(lnp);								    */
  /*     f+=lnp;								    */
  /*   }									    */
  /* }										    */
  /* f=-f;									    */
  /**********************************************************************************/

  for (i=0;i<N;++i){
    n1=n_i[i];
    n2=0.0;
    for (ii=0;ii<N;++ii){
      for (jj=0;jj<n_i[ii];++jj){
  	num=n_i[i]*exp(g_i[i])*exp_kii_xij[i][ii][jj];
  	din=0.0; for (iii=0;iii<N;++iii) din+=n_i[iii]*exp(g_i[iii])*exp_kii_xij[iii][ii][jj];
  
  	n2+=num/din;
      }
    }
    g[i]=-n1+n2;
  }

  /******************************************************************/
  /* for (i=0;i<N;++i) {					    */
  /*   g[i]=0.0;						    */
  /*   for (ii=0;ii<N;++ii) {					    */
  /*     for (jj=0;jj<n_i[ii];++jj) {				    */
  /* 	g[i]-=exp_kii_xij[i][ii][jj]*P_ij[ii][jj]*f_i[i];	    */
  /*     }							    */
  /*   }							    */
  /*   g[i]+=1.0;						    */
  /*   g[i]=n_i[i]*g[i];					    */
  /*   g[i]=-g[i];						    */
  /* }								    */
  /******************************************************************/

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

int optimize_lnL_BFGS(double *g, struct data dat ) {
  int i,j,k,ret = 0;
  lbfgsfloatval_t f;
  lbfgsfloatval_t *g_i = lbfgs_malloc(dat.K);
  lbfgs_parameter_t param;

  if (g_i == NULL) {
    printf("ERROR: Failed to allocate a memory block for variables.\n");
    return 1;
  }

  for (i=0;i<dat.K;++i) g_i[i] = log(1.0);

  //  for (i=0;i<dat.K;++i) w_ij[i] = log(w_ij[i]);
  
  lbfgs_parameter_init(&param);

  ret = lbfgs(dat.K, g_i, &f, evaluate, progress,/*(void *) */&dat, &param);

  printf("L-BFGS optimization terminated with status code = %d\n", ret);

  for (i=0;i<dat.K;++i) g[i] = g_i[i];

  lbfgs_free(g_i);

  return 0;
}
