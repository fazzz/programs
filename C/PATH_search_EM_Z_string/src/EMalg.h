
#ifndef INCLUDE_EMalg
#define INCLUDE_EMalg

int E_step(int N,int K,double **x_n,double **nyu_k,double ***Sigma_k,double *pi_k, double **gamma_nk,double pi);
double M_step(int N,int K,double **x_n,double **nyu_k,double ***Sigma_k,double *pi_k, double **gamma_nk);
double EM_lnrou(int N,int K,double **x_n,double **nyu_k,double ***Sigma_k,double *pi_k, double pi);

int E_step_fprob(int Nx, double minx,double dx,int Ny,double miny,double dy,
		 int K,double **prob,double **nyu_k,double ***Sigma_k,double *pi_k, double ***gamma_nk,
		 double pi);
double M_step_fprob(int Nx,double minx,double dx,int Ny, double miny,double dy,
		    /*int N,*/int K,double **prob,double **nyu_k,double ***Sigma_k,double *pi_k, double ***gamma_nk);
double EM_lnrou_fprob(int Nx, double minx,double dx, int Ny,double miny,double dy,
		      /*int N,*/int K,double **prob,double **nyu_k,double ***Sigma_k,double *pi_k, double pi);

double *Create_mixed_twoD_GaussianMap(double minx,double maxx,double dx,
				      double miny,double maxy,double dy,
				      double **nyu,double ***Sigma, double *pi_k, int K, double pi);

#endif
