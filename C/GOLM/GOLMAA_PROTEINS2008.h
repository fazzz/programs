
#ifndef INCLUDE_GOLMAA_P2008
#define INCLUDE_GOLMAA_P2008

#define UNIT 4.184070*100.0

#include "GOLMAA_PROTEINS2008_set.h"

double GOLMAA_PROTEINS2008_ff_calcff(double *crd, int numatom,struct potential_GOLMAA_PROTEINS2008 *ene);

double GOLMAA_PROTEINS2008_ff_calcff_wobaimp(double *crd, int numatom,struct potential_GOLMAA_PROTEINS2008 *ene);

double GOLMAA_PROTEINS2008_ff_calcff_debug(double *crd, int numatom,struct potential_GOLMAA_PROTEINS2008 *ene,int numstep);

double GOLMAA_PROTEINS2008_ff_calcBOND(double *crd,int numatom,double *p_b,double **f_b,double Kb,double *bon_equ,int **pairs,int numbond);
double GOLMAA_PROTEINS2008_ff_calcANGLE(double *crd,int numatom,double *p_a,double **f_a,double Ka,double *ang_equ,int **pairs,int numangl);
//double GOLMAA_PROTEINS2008_ff_calcDIHE(double *crd,int numatom,double *p_d,double **f_d,double Kd1,double Kd2,double *dih_equ,int **pairs,int numdihe);
//double GOLMAA_PROTEINS2008_ff_calcDIHE(double *crd,int numatom,double *p_d,double **f_d,double Kd1,double Kd2,double Kdi,double *dih_equ,int **pairs,int numdihe);
double GOLMAA_PROTEINS2008_ff_calcDIHE(double *crd,int numatom,double *p_d,double **f_d,double Kd1,double Kd2,double Kdi,double *dih_equ,int **pairs,int numdihe, int *impindex);

double GOLMAA_PROTEINS2008_ff_calcDIHE_debug(double *crd,int numatom,double *p_d,double **f_d,double Kd1,double Kd2,double Kdi,double *dih_equ,int **pairs,int numdihe, int *impindex,int numstep );

double GOLMAA_PROTEINS2008_ff_calcff_nonlocal(double *crd, int numatom,int *NC_index,int numNC,int *NotNC_index,int numNotNC,
					      double *ALJ_natatt,double *BLJ_natatt,double e_natatt,double ALJ_Repul,
					      double *p_natatt,double *p_repul,double **f_natatt,double **f_repul);

double GOLMAA_PROTEINS2008_ff_calcDIHE_woimp(double *crd,int numatom,double *p_d,double **f_d,double Kd1,double Kd2,double Kdi,double *dih_equ,int **pairs,int numdihe, int *impindex);

double GOLMAA_PROTEINS2008_ff_calcff_wobaimp_b(double *crd, int numatom,struct potential_GOLMAA_PROTEINS2008 *ene);

double GOLMAA_PROTEINS2008_ff_calcff_b(double *crd, int numatom,struct potential_GOLMAA_PROTEINS2008 *ene);

double GOLMAA_PROTEINS2008_ff_calcff_nonlocal_b(double *crd, int numatom,int *NC_index,int numNC,int *NotNC_index,int numNotNC,
						double *ALJ_natatt,double *BLJ_natatt,double e_natatt,double ALJ_Repul,
						double *p_natatt_t,double *p_repul_t,double **f_natatt,double **f_repul);

#endif
