
#ifndef INCLUDE_MD
#define INCLUDE_MD

#include <netcdf.h>
#include "netcdf_mineL.h"

#include "FFLc.h"

double runMDb_NHC_MP1998_Amber_AAFF(double *crd,double *vel, double *mass, int numatom,
				    double *zeta,double *V_zeta, double Q,
				    struct potential e, struct force f,
				    double T, double NfKT, int numstep, int interval,int *l,
				    double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
				    double *avePE, double *aveKE,double *aveT,
				    double *varPE, double *varKE,double *varT, double UNITT, double k_B,
				    struct my_netcdf_out_id_MCD nc_id_MCD,  FILE *outputfile,struct AmberParmL ap);

#endif
