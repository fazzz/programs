
#ifndef INCLUDE_GOLM_Clementi_MB
#define INCLUDE_GOLM_Clementi_MB

#include "GOLM_Clementi_set.h"

struct potential_GOLM_Clementi_MB {
  double p_MB;

  double **f_MB;

  struct potential_GOLM_Clementi e1;
  struct potential_GOLM_Clementi e2;

};

double GOLM_Clementi_MB_ff_calcff(double *crd, int numatom, double de, double d2,
				  struct potential_GOLM_Clementi_MB *ene);

double GOLM_Clementi_MB_ff_set_calcff(struct potential_GOLM_Clementi_MB *ene,
				      double *refcrd1,double *refcrd2,
				      double *refcrdAA1,double *refcrdAA2,
				      int numCAatom, int numatom, double ep);

double GOLM_Clementi_MB_Kai(double *crd,int numatom, double de,double d, double d2, 
			    struct potential_GOLM_Clementi_MB *ene);

double MD_Propagetor_NH_MP1998_GOLM_Clementi_MB(double *crd,double *vel,double *mass,
						double *zeta,double *V_zeta,double Q,double NfKT,
						int numatom,double *KEv,double *PEv,
						double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
						struct potential_GOLM_Clementi_MB *e_GOLM, double de, double d2);

#endif
