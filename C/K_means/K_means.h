
#ifndef INCLUDE_Kmeans
#define INCLUDE_Kmeans

double Kmeans_Estep(int N, int K, double **x, double **nyu, double **gamma_nk, double **dist);

double Kmeans_Mstep(int N, int K, double **x, double **nyu, double **gamma_nk);

double Kmeans_J(int N, int K,double **x,double **nyu, double **gamma_nk, double **dist);

double Kmeans_Estep_fprob(int Nx, double minx,double dx,int Ny,double miny,double dy,
			  int K, double **prob, double **nyu, double ***gamma_nk, double ***dist);

double Kmeans_Mstep_fprob(int Nx, double minx,double dx,int Ny,double miny,double dy,
			  int K, double **prob, double **nyu, double ***gamma_nk);

double Kmeans_J_fprob(int Nx, double minx,double dx,int Ny,double miny,double dy,
		      int K,double **prob,double **nyu, double ***gamma_nk, double ***dist);

#endif


