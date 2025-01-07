
#ifndef INCLUDE_GOLM_set
#define INCLUDE_GOLM_set

#define FC_bond_Clementi 100.0
#define FC_angl_Clementi 20.0
#define FC_dihe1_Clementi 10.0
#define FC_dihe2_Clementi 0.5

#define ep_natatt_Clementi 0.18
#define ep_repul_Clementi  0.18

#define cradii_repul_Clementi 4.0

int GOLMff_set_calcff(struct potential_GOLM *ene, double *cord,int numatom);

int *set_BOND_param(int *num_bond,double *refcrd);
int *set_ANGL_param(int *num_angl,double *refcrd);
int *set_DIHE_param(int *num_dihe,double *refcrd);
int *set_Repul_param(int *num_repul);

#endif
