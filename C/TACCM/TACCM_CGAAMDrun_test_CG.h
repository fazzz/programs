
#ifndef INCLUDE_TA_MD
#define INCLUDE_TA_MD

#include <netcdf.h>
#include "netcdf_mineL.h"

#include "FFL.h"

double ffL_calcffandforce_AA(double *crd, int numatom,struct potential *ene,struct force *f) ;

double TACCM_MD_Propagetor_NH_MP1998_CG_test(double *crd,double *vel,double *mass,
					     double *zeta,double *V_zeta,double Q,
					     double NfKT,int numatom,double *KE,double *KEv,double *PEv,
					     double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
					     double parameterCG, struct potential *e, struct force *f,
					     double *Z,  int numZ,double *theta,double Kapa,
					     int **pairs, double *PEZ,double pi);

double ffL_calcffandforce_CG(double *crd, int numatom,
			     struct potential *ene,struct force *f,double parameterCG);

int ffL_calcDIHE_CG(double *p_d,
		    double *n_d,
		    double *cord,
		    int flagp, int flagf, int flaginp,double parameterCG);

int ffL_calcDIHE_force_Cartesian_CG(double *f_d,double *cord,double parameterCG);

#endif
