
#ifndef INCLUDE_Gaussian
#define INCLUDE_Gaussian

double twoD_Gaussian(double *x,double *nyu,double **Sigma,double pi);

double de_twoD_Gaussian(double *x,double **f,double *nyu,double **Sigma,double pi);

double mixed_twod_Gaussian(double *x,double **nyu,double ***Sigma,double *pi_k,int K,double pi);

void de_ln_twoD_Mixed_Gaussian(double *x,double *f,double **nyu,double ***Sigma,double *pi_k,int K,double pi);

double *Create_twoD_GaussianMap(double minx,double maxx,double dx,
				double miny,double maxy,double dy,
				double *nyu,double **Sigma,double pi);

double *Create_mixed_twoD_GaussianMap(double minx,double maxx,double dx,
				      double miny,double maxy,double dy,
				      double **nyu,double ***Sigma,
				      double *pi_k,int K,double pi);

#endif

