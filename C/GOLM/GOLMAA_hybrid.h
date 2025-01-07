
#ifndef INCLUDE_GOLMAA_hb
#define INCLUDE_GOLMAA_hb

#define UNIT 4.184070*100.0

double GOLMAA_hyb_ff_calcff(double *crd, int numatom,struct potential_GOLMAA_hybrid *ene);

double GOLMAA_hyb_ff_calcff_nonlocal(double *crd, int numatom,int *NC_index,int numNC,int *NotNC_index,int numNotNC,
				     double *ALJ_natatt,double *BLJ_natatt,double e_natatt,double ALJ_Repul,
				     double *p_natatt,double *p_repul,double **f_natatt,double **f_repul);


#endif
