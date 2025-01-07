
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLM_Clementi_set_wmCutOff.h"
#include "GOLM_Clementi.h"
#include "GOLM_Clementi_MB_wmCutOff.h"

double GOLM_Clementi_MB_ff_set_calcff_wmCutOff(struct potential_GOLM_Clementi_MB *ene,
					       double *refcrd1,double *refcrd2,
					       double *refcrdAA1,double *refcrdAA2,
					       int numCAatom, int numatom, double ep, double cutoff) {
  int i,j;

  GOLM_Clementi_ff_set_calcff2_wmCutOff(&((*ene).e1),refcrd1,refcrdAA1,numCAatom,numatom,ep,cutoff);
  GOLM_Clementi_ff_set_calcff2_wmCutOff(&((*ene).e2),refcrd2,refcrdAA2,numCAatom,numatom,ep,cutoff);

  (*ene).f_MB = (double **)gcemalloc(sizeof(double *)*numCAatom);
  for (i=0;i<numCAatom;++i) (*ene).f_MB[i]=(double *)gcemalloc(sizeof(double)*3);

  return 0.0;
}
