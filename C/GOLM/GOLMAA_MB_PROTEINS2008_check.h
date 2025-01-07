

#ifndef INCLUDE_GOLMAA_MB_C_P2008
#define INCLUDE_GOLMAA_MB_C_P2008

double GOLMAA_MB_PROTEINS2008_calcff_check(char *inputfilename,
					   char *reffilename1,char *reffilename2,char *parmfilename,
					   int numspatom,double dx,
					   double d,double de,
					   double ep,int nibnum,double criteria,
					   double f_MB[3],double f_e1[3],double f_e2[3]);

#endif
