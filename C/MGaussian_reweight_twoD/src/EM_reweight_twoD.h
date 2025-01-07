
#ifndef INCLUDE_EM_reweight
#define INCLUDE_EM_reweight

void E_step_twoD_rw(int num_sim, int *n, double ***x,                     // # of simulation, # of snapshots, data,
		    double *minx, double *maxx, int num_Simpson,          // parameters for Simpson integration 1
		    double **k_umbrella, double **x0,                     // parameters for Simpson integration 2
		    int K, double **nyu_k,double ***Sigma_k,double *pi_k, // parameters of MG
		    double ***gammak_ij,  double **z_ik, double *z_i,     // responsibilities
		    double pi,int periodicflag, double periodicity);

void M_step_twoD_rw(int num_sim, int *n,int n_sum, double ***x,           // # of simulation, # of snapshots, data,
		    double *minx, double *maxx, int num_Simpson,          // parameters for Simpson integration 1
		    double **k_umbrella, double **x0,                     // parameters for Simpson integration 2
		    int K, double **nyu_k,double ***Sigma_k,double *pi_k, // parameters of MG
		    double ***gammak_ij,  double **z_ik, double *z_i,     // responsibilities
		    double pi, int periodicflag, double periodicity);

double EM_L_twoD_rw(int num_sim, int *n, double ***x,             // # of simulation, # of snapshots, data,
		    double *z_i,                                  // free energy differences bet. simulations
		    int K, double **nyu_k,double ***Sigma_k,double *pi_k, // parameters of MG
		    double pi);

#endif
