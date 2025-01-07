
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "lbfgs.h"
#include "Optimize_BFGS_GBRW_multiD_2.h"
#include "EF.h"

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *g_k,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step)
{
  int i,j,k,l,m,ii,jj,kk,ll;

  lbfgsfloatval_t f = 0.0;

  double n1,n2,n11,n22,n111,n222,n1111,n2222,num,din;

  double pi,sqrt_pi;

  int K,k_each,N,*n_i,dim;
  double *h_l,**x_l_i,***x_l_ij;
  double ***C_l_ik,****S_l_k_ij;
  double **A_l_k;

  //  double **g_l_k;
  double *g_K;

  struct data *dat = instance;

  K=(*dat).K;
  k_each=(*dat).k_each;
  N=(*dat).n_sim;
  dim=(*dat).dim;

  h_l=(double *)gcemalloc(sizeof(double)*(*dat).dim);
  for (i=0;i<(*dat).dim;++i) h_l[i]=(*dat).h_l[i];

  n_i=(int *)gcemalloc(sizeof(int)*N);
  for (i=0;i<N;++i) n_i[i]=(*dat).n[i];

  x_l_ij=(double ***)gcemalloc(sizeof(double **)*dim);
  for (l=0;l<(*dat).dim;++l) 
    x_l_ij[l]=(double **)gcemalloc(sizeof(double *)*N);
  for (l=0;l<(*dat).dim;++l) 
    for (i=0;i<N;++i) 
      x_l_ij[l][i]=(double *)gcemalloc(sizeof(double)*n_i[i]);
  for (l=0;l<(*dat).dim;++l) 
    for (i=0;i<N;++i) 
      for (j=0;j<n_i[i];++j) 
	x_l_ij[l][i][j]=(*dat).x_l_ij[l][i][j];

  C_l_ik=(double ***)gcemalloc(sizeof(double **)*dim);
  for (l=0;l<dim;++l)
    C_l_ik[l]=(double **)gcemalloc(sizeof(double *)*N);
  for (l=0;l<dim;++l)
    for (i=0;i<N;++i) 
      C_l_ik[l][i]=(double *)gcemalloc(sizeof(double *)*K);
  for (l=0;l<dim;++l)  
    for (i=0;i<N;++i) 
      for (k=0;k<K;++k) C_l_ik[l][i][k]=(*dat).C_l_ik[l][i][k];

  S_l_k_ij=(double ****)gcemalloc(sizeof(double ***)*dim);
  for (l=0;l<dim;++l)
    S_l_k_ij[l]=(double ***)gcemalloc(sizeof(double **)*K);
  for (l=0;l<dim;++l)
    for (i=0;i<K;++i) 
      S_l_k_ij[l][i]=(double **)gcemalloc(sizeof(double *)*N);
  for (l=0;l<dim;++l)
    for (i=0;i<K;++i)
      for (j=0;j<N;++j) 
	S_l_k_ij[l][i][j]=(double *)gcemalloc(sizeof(double)*n_i[j]);
  for (l=0;l<dim;++l)
    for (k=0;k<K;++k)
      for (i=0;i<N;++i)
	for (j=0;j<n_i[i];++j)
	  S_l_k_ij[l][k][i][j]=(*dat).S_l_k_ij[l][k][i][j];

  A_l_k=(double **)gcemalloc(sizeof(double *)*dim);
  for (l=0;l<dim;++l)
    A_l_k[l]=(double *)gcemalloc(sizeof(double)*K);
  for (l=0;l<dim;++l)
    for (k=0;k<K;++k)
      A_l_k[l][k]=(*dat).A_l_k[l][k];

  /*****************************************************/
  /* g_l_k=(double **)gcemalloc(sizeof(double *)*dim); */
  /* for (l=0;l<dim;++l)			       */
  /*   g_l_k[l]=(double *)gcemalloc(sizeof(double)*K); */
  /*****************************************************/

  g_K=(double *)gcemalloc(sizeof(double)*K);
  
  /***************************/
  /* i=0;		     */
  /* for (l=0;l<dim;++l) {   */
  /*   for (k=0;k<K;++k) {   */
  /*     g_l_k[l][k]=g_k[i]; */
  /*     ++i;		     */
  /*   }		     */
  /* }			     */
  /***************************/

  for (i=0;i<K;++i) {
    g_K[i]=g_k[i];
  }
    
  pi=acos(-1.0);

  f=0.0;
  n1=0.0;
  n2=0.0;
  for (i=0;i<N;++i){
    n11=0.0; 
    for (k=0;k<K;++k) {
      n111=1.0;
      for (l=0;l<dim;++l) {
	n111=n111*A_l_k[l][k]*C_l_ik[l][i][k]; 
      }
      n11+=exp(g_K[k])*n111;
    }
    n11=log(n11);
    n1+=n_i[i]*n11;
  }

  for (i=0;i<N;++i){
    for (j=0;j<n_i[i];++j){
      n22=0.0; 
      for (k=0;k<K;++k) {
	n222=1.0;
	for (l=0;l<dim;++l) {
	  n222=n222*A_l_k[l][k]*S_l_k_ij[l][k][i][j]; 
	}
	n22+=exp(g_K[k])*n222;
      }
      n22=log(n22);
      n2+=n22;
    }
  }
  f=n1-n2;

  for (k=0;k<K;++k){
    n1=0.0;
    for (i=0;i<N;++i) {
      n11=1.0;
      for (l=0;l<dim;++l) {
	n11=n11*A_l_k[l][k]*C_l_ik[l][i][k];
      }
      num=exp(g_K[k])*n11;

      din=0.0; 
      for (kk=0;kk<K;++kk) {
	n111=1.0;
	for (ll=0;ll<dim;++ll) {
	  n111=n111*A_l_k[ll][kk]*C_l_ik[ll][i][kk];
	}
	din+=exp(g_K[kk])*n111;
      }
      n1+=n_i[i]*num/din;
    }

    n2=0.0;
    for (i=0;i<N;++i){
      for (j=0;j<n_i[i];++j){
	n222=1.0;
	for (l=0;l<dim;++l) {
	  n222=n222*A_l_k[l][k]*S_l_k_ij[l][k][i][j];
	}
	num=exp(g_K[k])*n222;

	din=0.0; 
	for (kk=0;kk<K;++kk) {
	  n2222=1.0;
	  for (ll=0;ll<dim;++ll) {
	    n2222=n2222*A_l_k[ll][kk]*S_l_k_ij[ll][kk][i][j];
	  }
	  din+=exp(g_K[kk])*n2222;
	}
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
  
//int optimize_lnL_BFGS_multiD_2(double **g_l_k, struct data dat ) {
int optimize_lnL_BFGS_multiD_2(double *g, struct data dat ) {
  int i,j,k,ret = 0;
  lbfgsfloatval_t f;
  //  lbfgsfloatval_t *g_k = lbfgs_malloc(dat.K*dat.dim);
  lbfgsfloatval_t *g_k = lbfgs_malloc(dat.K);
  lbfgs_parameter_t param;

  if (g_k == NULL) {
    printf("ERROR: Failed to allocate a memory block for variables.\n");
    return 1;
  }

  /********************************/
  /* k=0;			  */
  /* for (i=0;i<dat.dim;++i) {	  */
  /*   for (j=0;j<dat.K;++j) {	  */
  /*     g_k[k] = log(1.0/dat.K); */
  /*     ++k;			  */
  /*   }			  */
  /* }				  */
  /********************************/

  for (i=0;i<dat.K;++i) {
    g_k[i] = log(1.0/dat.K);
  }

  lbfgs_parameter_init(&param);

  //  ret = lbfgs(dat.K*dat.dim, g_k, &f, evaluate, progress,&dat, &param);
  ret = lbfgs(dat.K, g_k, &f, evaluate, progress,&dat, &param);

  printf("L-BFGS optimization terminated with status code = %d\n", ret);

  /*****************************/
  /* k=0;		       */
  /* for (i=0;i<dat.dim;++i) { */
  /*   for (j=0;j<dat.K;++j) { */
  /*     g_l_k[i][j] = g_k[k]; */
  /*     ++k;		       */
  /*   }		       */
  /* }			       */
  /*****************************/

  for (i=0;i<dat.K;++i) {
    g[i] = g_k[i];
  }

  lbfgs_free(g_k);

  return 0;
}
