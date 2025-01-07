
#ifndef INCLUDE_GOLM_Cl_set
#define INCLUDE_GOLM_Cl_set

#define criteria_Clementi 6.5
#define ep_natatt_Clementi 0.36
#define epsilon_Clementi 1.0
#define cradii_repul_Clementi 4.0

struct potential_GOLM_Clementi {
  // energy terms
  double p_t;
  double p_natatt_t,p_repul_t,p_d_t,p_a_t,p_b_t;
  double *p_natatt,*p_repul,*p_d,*p_a,*p_b;

  // force terms
  double **f_t;
  double **f_natatt,**f_repul,**f_d,**f_a,**f_b;

  // parameters
  double *ALJ_natatt,*BLJ_natatt;
  double ep_natatt;
  double *cradii_natatt;
  int  num_natatt,*index_natatt;

  double ALJ_repul;
  double ep_repul;
  double cradii_repul;
  int num_repul,*index_repul;

  int **bpairs,**apairs,**dpairs;
  double *bon_equ,*ang_equ,*dih_equ;
  double Kb,Ka,Kd1,Kd2;
};

int GOLM_Clementi_ff_set_calcff(struct potential_GOLM_Clementi *ene, double *refcrd,double *refcrdAA, int numatom,int numatomAA);

int GOLM_Clementi_ff_set_calcff2(struct potential_GOLM_Clementi *ene, double *refcrd,double *refcrdAA, int numatom,int numatomAA, double ep);

int **GOLM_Clementi_make_native_contact(double *refcrdAA,double criteria,int *numnc,int numatom, int numCAatom);

#endif


