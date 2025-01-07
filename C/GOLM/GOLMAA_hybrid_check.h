
#ifndef INCLUDE_GOLMAA_hb_c
#define INCLUDE_GOLMAA_hb_c

#define UNIT 4.184070*100.0

double GOLMAAff_hybrid_calcff_check(double *crd, int numatom,
				    struct potential_GOLMAA_hybrid *ene,
				    int numspatom,double dx,
				    double f_natatt[3],double f_repul[3]);

#endif

