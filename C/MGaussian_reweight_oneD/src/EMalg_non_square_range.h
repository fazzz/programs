
#ifndef INCLUDE_EMalg_nsr
#define INCLUDE_EMalg_bsr

int E_step_fprob_nsr(int Nx, double minx,double dx,int Ny,double miny,double dy,
		     int K,double **prob,double **nyu_k,double ***Sigma_k,double *pi_k, double ***gamma_nk,
		     double pi);

double M_step_fprob_nsr(int Nx,double minx,double dx,int Ny, double miny,double dy,
			/*int N,*/int K,double **prob,double **nyu_k,double ***Sigma_k,double *pi_k, double ***gamma_nk);

double EM_lnrou_fprob_nsr(int Nx, double minx,double dx, int Ny,double miny,double dy,
			  /*int N,*/int K,double **prob,double **nyu_k,double ***Sigma_k,double *pi_k, double pi);

double EM_lnrou_fprob_nsr2(int Nx, double minx,double dx, int Ny,double miny,double dy,
			   /*int N,*/int K,double **prob,double **nyu_k,double ***Sigma_k,double *pi_k, double pi);

double M_step_fprob_nsr_2013_06_07(int Nx,double minx,double dx,int Ny, double miny,double dy,
				   /*int N,*/int K,double **prob,double **nyu_k,double ***Sigma_k,double *pi_k, double ***gamma_nk);


#endif
