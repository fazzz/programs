
#ifndef INCLUDE_GOLMAA_hb_set
#define INCLUDE_GOLMAA_hb_set

#define criteria_hybrid 5.0
#define ep_natatt_hybrid 3.8
#define ep_repul_hybrid 1.0
#define cradii_repul_hybrid 2.0

struct potential_GOLMAA_hybrid {
  // energy terms
  double p_t;

  double p_natatt_t,p_repul_t;
  double *p_natatt,*p_repul;

  // force terms
  double **f_t;
  double **f_natatt,**f_repul;

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

};

int GOLMAA_hybrid_ff_set_calcff(struct potential_GOLMAA_hybrid *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond);

int GOLMAA_hybrid_ff_set_calcff_wtune(struct potential_GOLMAA_hybrid *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond, double ep);

int **GOLMAA_hybrid_ff_set_make_native_contact(double *refcrd,double criteria,int *numNC,int numatom,int numres);


int GOLMAA_hybrid_ff_set_calcff_2(struct potential_GOLMAA_hybrid *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond);

int GOLMAA_hybrid_ff_set_calcff_2_wtune(struct potential_GOLMAA_hybrid *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond, double ep);

int **GOLMAA_hybrid_ff_set_make_native_contact_2(double *refcrd,double criteria,int *numNC,int numatom,int numres);

int GOLMAA_hybrid_ff_set_calcff_3(struct potential_GOLMAA_hybrid *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond);

int GOLMAA_hybrid_ff_set_calcff_3_wtune(struct potential_GOLMAA_hybrid *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond, double ep);

int **GOLMAA_hybrid_ff_set_make_native_contact_3(double *refcrd,double criteria,int *numNC,int numatom,int numres);

int GOLMAA_hybrid_ff_set_calcff_4_wtune(struct potential_GOLMAA_hybrid *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond, double ep);

int **GOLMAA_hybrid_ff_set_make_native_contact_4(double *refcrd,double criteria,int *numNC,int numatom,int numres);

int **GOLMAA_hybrid_ff_set_make_native_contact_5(double *refcrd,double criteria,int *numNC,int numatom,int numres,int nibnum);

int GOLMAA_hybrid_ff_set_calcff_6_wtune(struct potential_GOLMAA_hybrid *ene, double *refcrd, int numatom,int numres,int *non_bonding_index, int numnonbond, double ep, int nibnum,double criteria);

int **GOLMAA_hybrid_ff_set_make_native_contact_6(double *refcrd,double criteria,int *numNC,int numatom,int numres,int nibnum);

#endif


