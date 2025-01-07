#ifndef INCLUDE_FF_c
#define INCLUDE_FF_c

double ffL_calcffandforce_check(char *inputfilename,char *parmfilename,
				int numspatom,double dx,
				double f_es[3],double f_LJ[3],
				double f_14_es[3],double f_14_LJ[3],
				double f_d[3], double f_b[3], double f_a[3]
				);

double UMB_calcffandforce_check(double *crd,int numatom,
				int *pairp,int num,double *fcp,double *dih_equ,
				int numspatom, double dx, double f_UMB[3]);

#endif
