
#ifndef INCLUDE_GOLM_Clementi_MB_wmCO
#define INCLUDE_GOLM_Clementi_MB_wmCO

#include "GOLM_Clementi_set.h"
#include "GOLM_Clementi_MB.h"

double GOLM_Clementi_MB_ff_set_calcff_wmCutOff(struct potential_GOLM_Clementi_MB *ene,
					       double *refcrd1,double *refcrd2,
					       double *refcrdAA1,double *refcrdAA2,
					       int numCAatom, int numatom, double ep, double cutoff);

#endif
