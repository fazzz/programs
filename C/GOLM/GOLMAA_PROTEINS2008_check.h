
#ifndef INCLUDE_GOLMAA_P2008_check_set
#define INCLUDE_GOLMAA_P2008_check_set

double GOLMAA_PROTEINS2008_calcff_check(char *inputfilename,char *reffilename,char *parmfilename,
					int numspatom,double dx,
					double ep,int nibnum,double criteria,
					double f_natatt[3],double f_repul[3],
					double f_d[3],double f_a[3],double f_b[3],
					double f_d1[4][3],double f_d2[4][3], int nums,
					double f_as[3][3], int numas);

double GOLMAA_PROTEINS2008_ff_calcANGLE_forcheck(double *crd,int numatom,
						 double Ka,double *ang_equ,
						 int **pairs,int numangl,
						 double *p_a,double f_a[3][3],int nums);

#endif

