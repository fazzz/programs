
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lbfgs.h"
#include "Optimize_BFGS_TPBRW.h"
#include "Simpson_integ_TPBRW.h"
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
  double *x_i,**x_ij;
  double *C_i,/****/C_sin_ik,/****/C_cos_ik;
  double ***sin_k_ij,***cos_k_ij;
  double *A;

  struct data *dat = instance;

  int num_Simpson;
  double minv,maxv;
  double *bk_i,beta;
  double a0;
  double *a_k;
  double *b_k;

  K=(*dat).K;
  N=(*dat).n_sim;

  num_Simpson=(*dat).num_Simpson;
  minv=(*dat).minv;
  maxv=(*dat).maxv;

  a0=g_k[0];
  a_k=(double *)gcemalloc(sizeof(double)*K);
  for (i=0;i<K;++i) a_k[i]=g_k[i+1];
  b_k=(double *)gcemalloc(sizeof(double)*K);
  for (i=0;i<K;++i) b_k[i]=g_k[i+K+1];

  bk_i=(double *)gcemalloc(sizeof(double)*N);
  for (i=0;i<N;++i) bk_i[i]=(*dat).bk_i[i];
  beta=(*dat).beta;
  
  n_i=(int *)gcemalloc(sizeof(int)*N);
  for (i=0;i<N;++i) n_i[i]=(*dat).n[i];

  x_i=(double *)gcemalloc(sizeof(double)*N);
  for (i=0;i<N;++i) x_i[i]=(*dat).x_i[i];

  x_ij=(double **)gcemalloc(sizeof(double *)*N);
  for (i=0;i<N;++i) x_ij[i]=(double *)gcemalloc(sizeof(double)*n_i[i]);
  for (i=0;i<N;++i) for (j=0;j<n_i[i];++j) x_ij[i][j]=(*dat).x_ij[i][j];

  C_i=(double *)gcemalloc(sizeof(double)*N);
  //  for (i=0;i<N;++i) C_i[i]=(*dat).C_i[i];

  //  C_sin_ik=(double **)gcemalloc(sizeof(double *)*N);
  //  for (i=0;i<N;++i) C_sin_ik[i]=(double *)gcemalloc(sizeof(double *)*K);
  //  for (i=0;i<N;++i) for (j=0;j<K;++j) C_sin_ik[i][j]=(*dat).C_sin_ik[i][j];

  //  C_cos_ik=(double **)gcemalloc(sizeof(double *)*N);
  //  for (i=0;i<N;++i) C_cos_ik[i]=(double *)gcemalloc(sizeof(double *)*K);
  //  for (i=0;i<N;++i) for (j=0;j<K;++j) C_cos_ik[i][j]=(*dat).C_cos_ik[i][j];

  sin_k_ij=(double ***)gcemalloc(sizeof(double **)*K);
  for (i=0;i<K;++i) sin_k_ij[i]=(double **)gcemalloc(sizeof(double *)*N);
  for (i=0;i<K;++i) for (j=0;j<N;++j) sin_k_ij[i][j]=(double *)gcemalloc(sizeof(double)*n_i[j]);
  for (k=0;k<K;++k)
    for (i=0;i<N;++i)
      for (j=0;j<n_i[i];++j)
	sin_k_ij[k][i][j]=(*dat).sin_k_ij[k][i][j];

  cos_k_ij=(double ***)gcemalloc(sizeof(double **)*K);
  for (i=0;i<K;++i) cos_k_ij[i]=(double **)gcemalloc(sizeof(double *)*N);
  for (i=0;i<K;++i) for (j=0;j<N;++j) cos_k_ij[i][j]=(double *)gcemalloc(sizeof(double)*n_i[j]);
  for (k=0;k<K;++k)
    for (i=0;i<N;++i)
      for (j=0;j<n_i[i];++j)
	cos_k_ij[k][i][j]=(*dat).cos_k_ij[k][i][j];

  pi=acos(-1.0);

  f=0.0;
  n1=0.0;
  n2=0.5*a0;
  for (i=0;i<N;++i){
    C_i[i]=Simpson_integ_oneD_TPBRW_C(num_Simpson,minv,maxv,x_i[i],K,bk_i[i],beta,a0,a_k,b_k);
    n1+=n_i[i]*log(C_i[i]);

    for (j=0;j<n_i[i];++j){
      for (k=0;k<K;++k)
	n2+=a_k[k]*sin_k_ij[k][i][j]+b_k[k]*cos_k_ij[k][i][j];
    }
  }
  n2=beta*n2;
  f=n1+n2;

  //  g[0]=0.0;
  n1=0.0;
  n2=0.0;
  for (i=0;i<N;++i){
    n1+=n_i[i]*0.5*beta;
    for (j=0;j<n_i[i];++j){
      n2+=0.5*beta;
    }
  }
  g[0]=-n1+n2;

  for (k=0;k<K;++k){
    n1=0.0;
    for (i=0;i<N;++i){
      C_sin_ik=Simpson_integ_oneD_TPBRW_C_sin(num_Simpson,minv,maxv,x_i[i],k,K,bk_i[i],beta,a0,a_k,b_k);
      n1+=n_i[i]*C_sin_ik/C_i[i];
    }

    n2=0.0;
    for (i=0;i<N;++i){
      for (j=0;j<n_i[i];++j){
	n2+=sin_k_ij[k][i][j];
      }
    }
    g[k+1]=beta*(-n1+n2);
  }

  for (k=0;k<K;++k){
    n1=0.0;
    for (i=0;i<N;++i){
      C_cos_ik=Simpson_integ_oneD_TPBRW_C_cos(num_Simpson,minv,maxv,x_i[i],k,K,bk_i[i],beta,a0,a_k,b_k);
      n1+=n_i[i]*C_cos_ik/C_i[i];
    }

    n2=0.0;
    for (i=0;i<N;++i){
      for (j=0;j<n_i[i];++j){
	n2+=cos_k_ij[k][i][j];
      }
    }
    g[k+K+1]=beta*(-n1+n2);
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

int optimize_lnL_BFGS_TPBRW(double *g, struct data dat ) {
  int i,j,k,ret = 0;
  lbfgsfloatval_t f;
  lbfgsfloatval_t *g_k = lbfgs_malloc(dat.K*2+1);
  lbfgs_parameter_t param;

  int N_bin=400;       // test
  double x,free_ene,p; // test
  double minx, maxx;   // test
  FILE *pmffile;       // test

  minx=dat.minv; // test
  maxx=dat.maxv; // test

  if (g_k == NULL) {
    printf("ERROR: Failed to allocate a memory block for variables.\n");
    return 1;
  }

  g_k[0] = 0.0;
  for (i=0;i<dat.K;++i) g_k[i+1]   = 1.0/dat.K;
  for (i=0;i<dat.K;++i) g_k[i+dat.K+1] = 1.0/dat.K;
  
  lbfgs_parameter_init(&param);

  ret = lbfgs((dat.K*2+1), g_k, &f, evaluate, progress,&dat, &param);

  // test
  pmffile=efopen("pmf.txt","w");
  for (i=0;i<N_bin;++i) {
    x=(maxx-minx)/(double)N_bin*i+minx;
    free_ene=0.5*g_k[0];
    for (j=0;j<dat.K;++j) {
      free_ene+=g_k[j+1]*sin((j+1)*x)+g_k[j+dat.K+1]*cos((j+1)*x);
    }
    fprintf(pmffile,"%10.8lf %10.8lf \n",x,free_ene);
  }
  fclose(pmffile);
  // test

  printf("L-BFGS optimization terminated with status code = %d\n", ret);

  for (i=0;i<dat.K*2+1;++i) g[i] = g_k[i];

  lbfgs_free(g_k);

  return 0;
}

