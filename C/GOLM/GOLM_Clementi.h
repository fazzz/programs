
#ifndef INCLUDE_GOLMAA_Cl
#define INCLUDE_GOLMAA_Cl

#include "GOLM_Clementi_set.h"

#define UNIT 4.184070*100.0

double GOLM_Clementi_ff_calcff(double *crd, int numatom,struct potential_GOLM_Clementi *ene);

double GOLM_Clementi_ff_calcff_natatt(double *crd, int numatom,int *index_numatt,int numatt,
				      double *ALJ_natatt,double *BLJ_natatt,double ep_natatt,
				      double *p_natatt,double **f_natatt);

double GOLM_Clementi_ff_calcff_non_natatt(double *crd, int numatom,
					  int *index_notnumatt,int numnotatt,
					  double ALJ_repul,
					  double *p_repul,double **f_repul);

double GOLM_Clementi_ff_calcBOND(double *crd,int numatom,double *p_b,double **f_b,double Kb,double *bon_equ);

double GOLM_Clementi_ff_calcANGLE(double *crd,int numatom,double *p_a,double **f_a,double Ka,double *ang_equ);

double GOLM_Clementi_ff_calcDIHE(double *crd,int numatom,double *p_d,double **f_d,double Kd1,double Kd2,double *dih_equ);

#endif
