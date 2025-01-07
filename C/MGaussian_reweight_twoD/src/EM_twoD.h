
#ifndef INCLUDE_EM
#define INCLUDE_EM

void E_step_twoD(int num_sim, int *n, double ***x,                     // # of simulation, # of snapshots, data,
		 int K, double **nyu_k,double ***Sigma_k,double *pi_k, // parameters of MG
		 double ***gammak_ij,                                  // responsibilities
		 double pi);

void M_step_twoD(int num_sim, int *n,int n_sum, double ***x,           // # of simulation, # of snapshots, data,
		 int K, double **nyu_k,double ***Sigma_k,double *pi_k, // parameters of MG
		 double ***gammak_ij,                                  // responsibilities
		 double pi);

double EM_L_twoD(int num_sim, int *n, double ***x,                      // # of simulation, # of snapshots, data,
		 int K, double **nyu_k,double ***Sigma_k,double *pi_k,  // parameters of MG
		 double pi);

#endif
