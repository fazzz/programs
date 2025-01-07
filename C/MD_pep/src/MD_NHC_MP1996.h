#ifndef INCLUDE_MD_MP1996
#define INCLUDE_MD_MP1996

#include "FFL.h"

#define ON 1
#define OFF 0

#define NVE 0
#define NVT 1

double MD_Generate_inivelo(double *vel,double *mass,int numatom,double KbT);

double MD_Propagetor_NH_Single_set_MP1996(int nc,double dt,double *dt2,double wdt2[3],double wdt4[3]);

double MD_Propagetor_NH_Single_part_MP1996(double *vel,double *mass,
					   double *zeta,double *V_zeta,double Q,
					   double NfKT,int numatom,int nc,
					   double wdt4[3],double wdt2[3]);

double MD_Propagetor_NH_MP1998_AAFF_Amber(double *crd,double *vel,double *mass,
					  double *zeta,double *V_zeta,double Q,
					  double NfKT,int numatom,double *KEv,double *PEv,
					  double dt,double dt2,int nc,double wdt4[3],double wdt2[3], 
					  struct potential *e, struct force *f);

double MD_Propagetor_vV_NVE_AAFF_Amber(double *crd,double *vel,double *mass,
				       int numatom,double dt,
				       struct potential *e, struct force *f);

#endif
