#ifndef INCLUDE_GBR_BFGS
#define INCLUDE_GBR_BFGS

#include "lbfgs.h"

struct data {
  double **x_ij;
  double  *x_i;
  double  *x_k;

  double *k;

  double beta;

  //  double *C_i,**C_sin_ik,**C_cos_ik;

  //  double *A;

  double ***sin_k_ij,***cos_k_ij;

  int *n;
  int n_sim;
  int K;

  int num_Simpson;
  double minv,maxv;

  double *bk_i;
};

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int K,
    const lbfgsfloatval_t step);

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
    int ls);

int optimize_lnL_BFGS_TPBRW(double *g, struct data dat );

#endif
