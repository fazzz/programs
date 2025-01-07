
#ifndef INCLUDE_Simpson_integ
#define INCLUDE_Simpson_integ

double Simpson_integ_C_oneD_Gaussian(int N, double minx,double maxx,
				     double k, double x0,double nyu,double Sigma,double pi);

double Simpson_integ_C_x_oneD_Gaussian(int N, double minx,double maxx,
				       double k, double x0,double nyu,double Sigma,double pi);

double Simpson_integ_C_x2_oneD_Gaussian(int N, double minx,double maxx,
					double k, double x0,double nyu,double Sigma,double pi);

double Simpson_integ_oneD_Gaussian(int N, double minx,double maxx,double nyu,double Sigma,double pi);

double oneD_Gaussian(double x,double nyu,double Sigma,double pi);

#endif
