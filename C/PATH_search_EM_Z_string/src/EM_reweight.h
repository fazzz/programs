
#ifndef INCLUDE_EM_reweight
#define INCLUDE_EM_reweight

void E_step_oneD(int num_sim, int *n, double **x,                      // # of simulation, # of snapshots, data,
		 double minx, double maxx, int num_Simpson,            // parameters for Simpson integration 1
		 double *k_umbrella, double *x0,                       // parameters for Simpson integration 2
		 int K, double *nyu_k,double *Sigma_k,double *pi_k,    // parameters of MG
		 double ***gammak_ij,  double **inv_f_ik, double *f_i, // responsibilities
		 double pi);

void M_step_oneD(int num_sim, int *n,int n_sum, double **x,            // # of simulation, # of snapshots, data,
		 double minx, double maxx, int num_Simpson,            // parameters for Simpson integration 1
		 double *k_umbrella, double *x0,                       // parameters for Simpson integration 2
		 int K, double *nyu_k,double *Sigma_k,double *pi_k,    // parameters of MG
		 double ***gammak_ij,  double **inv_f_ik, double *f_i, // responsibilities
		 double pi);

double EM_L(int num_sim, int *n, double **x,                    // # of simulation, # of snapshots, data,
	    double *f_i,                                        // free energy differences bet. simulations
	    int K, double *nyu_k,double *Sigma_k,double *pi_k,  // parameters of MG
	    double pi);

#endif
