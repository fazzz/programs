
#include <stdio.h>
#include <math.h>

#include "WHAM.h"
#include "EF.h"

static lbfgsfloatval_t evaluate(
    struct data *dat,
    const lbfgsfloatval_t *gs,
    lbfgsfloatval_t *grA,
    const int nv,
    const lbfgsfloatval_t step
    )
{
    int i,j,k,l,m,nn;
    double n1,n2,din,num;
    lbfgsfloatval_t A = 0.0;

    n1=0.0;
    n2=0.0;
    for (j = 0; j < nv+1; ++j) {
      for (nn = 0; nn < (*dat).n[j]; ++nn) {
	din=((*dat).n[0])*((*dat).expene[0][j][nn]);
	for (k=1;k<nv+1;++k) din+=((*dat).n[k])*((*dat).expene[k][j][nn])*exp(gs[k-1]);
	n2+=log(1.0/din);
      }
    }
    for (i=0;i<nv+1;++i)  n1+=((*dat).n[i])*gs[i-1];
    A=-n1-n2;

    for (i=1;i<nv+1;++i) {
      num=0.0;
      for (j=0;j<nv+1;++j) {
	for (nn = 0; nn < (*dat).n[j]; ++nn) {
	  din=(*dat).n[0]*((*dat).expene[0][j][nn]);
	  for (k=1;k<nv+1;++k) din+=((*dat).n[k])*((*dat).expene[k][j][nn])*exp(gs[k-1]);
	  num+=((*dat).expene[i][j][nn])/din;
	}
      }
      grA[i-1]=(*dat).n[i]*(exp(gs[i-1])*num-1.0);
    }

    return A;
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
  FILE *log;

  log=efopen("logWHAM.txt","w");
  fprintf(log,"Iteration %d:\n", k);
  fprintf(log,"  fx = %f\n", fx);
  for (i=0;i<n;++i)
    fprintf(log,"  x[%2d] = %f, \n",i, x[i]);
  fprintf(log,"  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
  fprintf(log,"\n");
  fclose(log);
  return 0;
}

int WHAM_fast_BFGS(double *expF, double ***expene, int n_sim, int *n) {
  int i,j,k,ret = 0;
  lbfgsfloatval_t A;
  lbfgsfloatval_t *g = lbfgs_malloc(n_sim-1);
  lbfgs_parameter_t param;

  struct data dat;

  dat.expene=(double ***)gcemalloc(sizeof(double **)*n_sim);
  for (i=0;i<n_sim;++i) dat.expene[i]=(double **)gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<n_sim;++i) for (j=0;j<n_sim;++j)	dat.expene[i][j]=(double *)gcemalloc(sizeof(double )*n[j]);

  dat.n=(int *)gcemalloc(sizeof(int)*n_sim);

  for (i=0;i<n_sim;++i) 
    for (j=0;j<n_sim;++j)	
      for (k=0;k<n[j];++k)	
	dat.expene[i][j][k]=expene[i][j][k];
  for (i=0;i<n_sim;++i) dat.n[i]=n[i];

  if (g == NULL) {
    printf("ERROR: Failed to allocate a memory block for variables.\n");
    return 1;
  }

  for (i=0;i<n_sim-1;++i) g[i] = 0.0;
  
  lbfgs_parameter_init(&param);

  ret = lbfgs(n_sim-1, g, &A, evaluate, progress, &dat, &param);

  //  printf("L-BFGS optimization terminated with status code = %d\n", ret);
  expF[0]=1.0;
  for (i=1;i<n_sim;++i) expF[i]=exp(g[i-1]);

  lbfgs_free(g);
  return 0;
}
