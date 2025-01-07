
#ifndef INCLUDE_GOLMAA_P2008_set
#define INCLUDE_GOLMAA_P2008_set

#define criteria_PROTEINS2008 6.5
#define ep_natatt_PROTEINS2008 0.36
#define epsilon_PROTEINS2008 1.0

#define ep_bond_PROTEINS2008 100
#define ep_angl_PROTEINS2008 20
#define ep_impd_PROTEINS2008 10
#define ep_dih1_PROTEINS2008 1.0
#define ep_dih2_PROTEINS2008 0.5

#define ep_repul_PROTEINS2008 0.01
#define cradii_repul_PROTEINS2008 2.5

struct potential_GOLMAA_PROTEINS2008 {
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

  int *NC_index,numNC;
  int *NotNC_index,numNotNC;

  double ALJ_repul;
  double ep_repul;
  double cradii_repul;
  int num_repul,*index_repul;

  int num_bond,num_angl,num_dihe;
  int **pairs_bond,**pairs_angl,**pairs_dihe;
  double *bon_equ,*ang_equ,*dih_equ;
  double Kb,Ka,Ki,Kd1,Kd2;
  int *impindex;
};

int GOLMAA_PROTEINS2008_ff_set_calcff(struct potential_GOLMAA_PROTEINS2008 *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond, double ep, int nibnum,double criteria);

int GOLMAA_PROTEINS2008_ff_set_calcff_b(struct potential_GOLMAA_PROTEINS2008 *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond, double ep, int nibnum,double criteria);

int **GOLMAA_PROTEINS2008_ff_set_make_native_contact(double *refcrd,double criteria,int *numNC,int numatom,int numres,int nibnum);

#endif


