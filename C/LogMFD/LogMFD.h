
#ifndef INCLUDE_LogMFD
#define INCLUDE_LogMFD

double LogMFD_Propagetor_NH_MP1998_1(double *crd, double *vel,double *mass,
				     double *zeta,double *V_zeta,double Q,
				     double NfKT,int N,double *KE,double *KEv,double *PEv,
				     double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
				     double *dF_dX, double F, double alpha, double gamma);

double LogMFD_Propagetor_NH_MP1998_2(double *crd, double *vel,double *mass,
				     double *zeta,double *V_zeta,double Q,
				     double NfKT,int N,double *KE,double *KEv,double *PEv,
				     double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
				     double *dF_dX, double F, double alpha, double gamma);

double LogMFD_Propagetor_NH_Single_part_MP1996(double *vel,double *mass,
					       double *zeta,double *V_zeta,double Q,
					       double NfKT,int N,int nc,double wdt4[3],double wdt2[3]);

double LogMFD_cF(double alpha, double gamma, double HMFD,
		 double NfKT, double KE,double KEv,double PEv);

#endif
