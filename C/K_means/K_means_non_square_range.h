
#ifndef INCLUDE_Kmeans_nsr
#define INCLUDE_Kmeans_nsr

double Kmeans_Estep_nsr(int N, int K, double **x, double **nyu, double **gamma_nk, double **dist);

double Kmeans_Mstep_nsr(int N, int K, double **x, double **nyu, double **gamma_nk);

double Kmeans_J_nsr(int N, int K,double **x,double **nyu, double **gamma_nk, double **dist);

double Kmeans_Estep_fprob_nsr(int Nx, double minx,double dx,int Ny,double miny,double dy,
			      int K, double **prob, double **nyu, double ***gamma_nk, double ***dist);

double Kmeans_Mstep_fprob_nsr(int Nx, double minx,double dx,int Ny,double miny,double dy,
			      int K, double **prob, double **nyu, double ***gamma_nk);

double Kmeans_J_fprob_nsr(int Nx, double minx,double dx,int Ny,double miny,double dy,
			  int K,double **prob,double **nyu, double ***gamma_nk, double ***dist);

#endif


