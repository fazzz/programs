#ifndef INCLUDE_WHAM
#define INCLUDE_WHAM

#include "lbfgs.h"

struct data {
  double ***expene;
  int *n;
};

static lbfgsfloatval_t evaluate(
    struct data *dat,
    const lbfgsfloatval_t *expf,
    lbfgsfloatval_t *grA,
    const int n_sim,
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

int WHAM_fast_BFGS(double *expF, double ***expene, int n_sim, int *n);

#endif
