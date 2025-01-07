

#ifndef INCLUDE_GOLMAA_P2008
#define INCLUDE_GOLMAA_P2008

#include "GOLMAA_PROTEINS2008_set.h"

double GOLMAA_PROTEINS2008_Amber_hybrid_ff_calcff_b(double *crd, int numatom,
						    struct potential_GOLMAA_PROTEINS2008 *ene);

double GOLMAA_PROTEINS2008_Amber_hybrid_calcDIHE_force_Cartesian(double **f_d,double *cord);
double GOLMAA_PROTEINS2008_Amber_hybrid_calcANGLE_force_Cartesian(double **f_a,double *cord);
double calcANGKE_force2(double atomi[3],double atomj[3],double atomk[3],double kang,double ang_eq,double *f);
double GOLMAA_PROTEINS2008_Amber_hybrid_calcBOND_force_Cartesian(double **f_b,double *cord);


#endif
