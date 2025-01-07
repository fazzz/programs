
#ifndef INCLUDE_Simpson_integ_SIGBRW
#define INCLUDE_Simpson_integ_SIGBRW


double Simpson_integ_oneD_TPBRW_C(int N, double minx,double maxx,double x_i,int K,double bk,double beta,
				  double a0, double *a_k, double *b_k);

double Simpson_integ_oneD_TPBRW_C_sin(int N, double minx,double maxx,double x_i,int k,int K,double bk,double beta,
				      double a0, double *a_k, double *b_k);

double Simpson_integ_oneD_TPBRW_C_cos(int N, double minx,double maxx,double x_i,int k,int K,double bk,double beta,
				      double a0, double *a_k, double *b_k);

#endif
