
#ifndef INCLUDE_GOLMAA_MB_P2008
#define INCLUDE_GOLMAA_MB_P2008

#include "GOLMAA_PROTEINS2008_set.h"

struct potential_GOLMAA_MB_PROTEINS2008 {
  double p_MB;

  double **f_MB;

  struct potential_GOLMAA_PROTEINS2008 e1;
  struct potential_GOLMAA_PROTEINS2008 e2;

};

double GOLMAA_MB_PROTEINS2008_ff_calcff(double *crd, int numatom,double de, double d2,
					struct potential_GOLMAA_MB_PROTEINS2008 *ene);

double GOLMAA_MB_PROTEINS2008_ff_calcff_set(struct potential_GOLMAA_MB_PROTEINS2008 *ene, double *refcrd1,double *refcrd2, int numatom,int numres,int *non_bonding_index, int numnonbond, double ep, int nibnum,double criteria);

double GOLMAA_MB_PROTEINS2008_Kai(double *crd,int numatom, double de,double d, double d2, struct potential_GOLMAA_MB_PROTEINS2008 *ene);

double GOLMAA_MB_PROTEINS2008_ff_calcff_wobaimp(double *crd, int numatom, double de, double d2,
						struct potential_GOLMAA_MB_PROTEINS2008 *ene);

double GOLMAA_MB_PROTEINS2008_wobaimp_Kai_debug(double *crd,int numatom, double de,double d, double d2, 
						struct potential_GOLMAA_MB_PROTEINS2008 *ene);

double GOLMAA_MB_PROTEINS2008_wobaimp_Kai(double *crd,int numatom, double de,double d, double d2, struct potential_GOLMAA_MB_PROTEINS2008 *ene);

double GOLMAA_MB_PROTEINS2008_ff_calcff_harmo(double x, 
					      double de, double d2,
					      double k1, double x1,
					      double k2, double x2,
					      double *f);

double GOLMAA_MB_PROTEINS2008_harmo_Kai(double x, 
					double de, double d,double d2,
					double k1, double x1,
					double k2, double x2);

#endif
