
#ifndef INCLUDE_AAFF_check
#define INCLUDE_AAFF_check

double AAFF_Amber_calcff_check(char *inputfilename,char *parmfilename,
			       int numspatom,double dx,
			       double f_e[3],double f_LJ[3],
			       double f_e_14[3],double f_LJ_14[3],
			       double f_d[3],double f_a[3],double f_b[3]);

#endif
