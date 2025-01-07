#ifndef INCLUDE_GBR_BFGS_multiD
#define INCLUDE_GBR_BFGS_multiD

#include "lbfgs.h"

struct data {
  double ***x_l_ij; // data of simulation ( # of dim. x # of simulation x # of # of step of each simulation )
  double **x_l_i;   // center of umbrella potential ( # of dim. x # of simulation )
  double **x_l_k;   // center of Gaussian base ( # of dim. x # of bases )

  double **k_l_i; // the force constant for each umbrella simulation ( # of dim. x # of simulation )
  double *h_l;    // the width of Gaussian of each axis ( # of dim. )

  double beta;

  double ***C_l_ik; // Int e^(xl-xli)^2*e^(xl-xkl)^2 dxl ( dim x n_sim x K )

  double **A_l_k; 

  double ****S_l_k_ij; // exp(x_l_ij-x_l_k) (dim x K x n_sim x n[i] )

  int *n;     // # of snapshots of each simulation
  int n_sim;  // # of simulations
  int dim;    // # of dimension
  int k_each; // # of bases along each axis
  int K;      // # of bases = dim x k_each
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

//int optimize_lnL_BFGS_multiD_2(double **g_l_k, struct data dat );
int optimize_lnL_BFGS_multiD_2(double *g_k, struct data dat );

#endif

