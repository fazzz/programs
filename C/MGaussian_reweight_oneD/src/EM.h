
#ifndef INCLUDE_EM
#define INCLUDE_EM

void E_step_oneD(int num_sim, int *n, double **x,                      // # of simulation, # of snapshots, data,
		 double minx, double maxx, int num_Simpson,            // parameters for Simpson integration 1
		 int K, double *nyu_k,double *Sigma_k,double *pi_k,    // parameters of MG
		 double ***gammak_ij,                                  // responsibilities
		 double pi,int periodicflag, double periodicity);

void M_step_oneD(int num_sim, int *n,int n_sum, double **x,            // # of simulation, # of snapshots, data,
		 double minx, double maxx, int num_Simpson,            // parameters for Simpson integration 1
		 int K, double *nyu_k,double *Sigma_k,double *pi_k,    // parameters of MG
		 double ***gammak_ij,                                  // responsibilities
		 double pi, int periodicflag, double periodicity);

double EM_L(int num_sim, int *n, double **x,                    // # of simulation, # of snapshots, data,
	    int K, double *nyu_k,double *Sigma_k,double *pi_k,  // parameters of MG
	    double pi);

#endif
