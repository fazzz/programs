
#ifndef INCLUDE_GOLM_check
#define INCLUDE_GOLM_check

double *GOLMff_calcff_check(double *crd, int numatom,struct potential_GOLM *ene,
			    int flagb,int flaga, int flagd, int flagnc, int flagnn,
			    int numspatom,double dx);

double *GOLMAAff_calcff_check(double *crd, int numatom,struct potential_GOLMAA *ene,
			      int numspatom,double dx,double f_natatt[3],double f_repul[3]);

#endif
