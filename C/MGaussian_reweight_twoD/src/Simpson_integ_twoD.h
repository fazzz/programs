
#ifndef INCLUDE_Simpson_integ
#define INCLUDE_Simpson_integ

#define ON 0
#define OFF 1

double Simpson_integ_twoD_Gaussian(int N, double *minx,double *maxx,double *nyu,double **Sigma,double pi);

double Simpson_integ_C_twoD_Gaussian(int N, double *minx,double *maxx,
				     double *k, double *x0,double *nyu,double **Sigma,double pi,
				     int periodicflag, double periodicity);

double *Simpson_integ_C_x_uk_twoD_Gaussian(int N, double *minx,double *maxx,
					   double *k, double *x0,double *nyu,double **Sigma,double pi,
					   int periodicflag, double periodicity);

double **Simpson_integ_C_x_uk2_twoD_Gaussian(int N, double *minx,double *maxx,
					     double *k, double *x0,double *nyu,double **Sigma,double pi,
					     int periodicflag, double periodicity);

#endif
