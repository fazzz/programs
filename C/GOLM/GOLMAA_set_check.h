
#ifndef INCLUDE_GOLMAA_set
#define INCLUDE_GOLMAA_set

#include "GOLMAA.h"

#define FC_dihe_Sanbommatsu 10.0

#define ep_natatt_Sanbonmatsu 1.0
#define ep_repul_Sanbonmatsu  0.01

#define cradii_repul_Sanbonmatsu 2.5

int GOLMAAff_set_calcff(struct potential_GOLMAA *ene, double *refcrd,int numatom,int **nb_matrix,double R_C_D,double constant);

//int *set_DIHE_param(int *num_dihe,double *refcrd);
//int *set_Repul_param(int *num_repul);

#endif
