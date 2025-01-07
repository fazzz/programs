
#ifndef INCLUDE_Simpson_integ_GB
#define INCLUDE_Simpson_integ_GB

double oneD_Gaussian_base(double x,double nyu,double Sigma);

double Simpson_integ_oneD_Gaussian_base(int N, double minx,double maxx,
					double nyu,double Sigma);

double Simpson_integ_C_oneD_Gaussian_base(int N, double minx,double maxx,
					  double k, double x0,double nyu,double Sigma,
					  int periodicflag, double periodicity);

double Simpson_integ_C_oneD_Gaussian_base(int N, double minx,double maxx,
					  double k, double x0,double nyu,double Sigma,
					  int periodicflag, double periodicity);

double Simpson_integ_C_P_GB(int N, double minx,double maxx,
			    int K, double *omega_k, 
			    double k, double x0, 
			    double *nyu,double *Sigma,
			    int periodicflag, double periodicity);

double C_func(double k, double x0, double x, int periodicflag, double periodicity);

#endif
