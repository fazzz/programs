#ifndef INCLUDE_MD_MP1996
#define INCLUDE_MD_MP1996

#include "GOLM_Clementi_set.h"
#include "GOLM_Clementi.h"
#include "GOLMAA_hybrid_set.h"
#include "GOLMAA_hybrid.h"
#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"
#include "GOLMAA_PROTEINS2008_debug.h"
#include "GOLMAA_MB_PROTEINS2008.h"
#include "UMBP.h"
#include "FFL.h"

double MD_Propagetor_NH_Single_set_MP1996(int nc,double dt,double *dt2,double wdt2[3],double wdt4[3]);
double MD_Propagetor_NH_Single_part_MP1996(double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,int nc,double wdt4[3],double wdt2[3]);
double MD_Propagetor_NH_MP1998_GOLM_Clementi(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLM_Clementi *e_GOLM);
double MD_Propagetor_NH_MP1998_GOLMAA_JCTC2011(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f, struct potential_GOLMAA_hybrid *e_GOLM);
double MD_Propagetor_NH_MP1998_14LJdab(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f);
double MD_Propagetor_NH_MP1998_GOLMAA_JCTC2011_wonat(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f, struct potential_GOLMAA_hybrid *e_GOLM);
double MD_Propagetor_NH_MP1998_14LJdab(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f);
double MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLMAA_PROTEINS2008 *e_GOLM);
double MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008_debug(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLMAA_PROTEINS2008 *e_GOLM,int bflag, int aflag, int dflag, int nflag);

double MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008_debug2(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLMAA_PROTEINS2008 *e_GOLM,int numstep);

double MD_Propagetor_NH_MP1998_AAFF_Amber(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f);
double MD_Propagetor_NH_MP1998_AAFF_Amber_wflag(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f,int bodflag,int angflag,int dihflag,int LJ14flag,int es14flag,int LJflag,int esflag);
double MD_Propagetor_NH_MP1998_AAFF_Amber_woH_wflag(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f,int bodflag,int angflag,int dihflag,int LJ14flag,int es14flag,int LJflag,int esflag);

double MD_Propagetor_NH_MP1998_AAFF_Amber_UMB_wflag(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f,int bodflag,int angflag,int dihflag,int LJ14flag,int es14flag,int LJflag,int esflag,int *pairp,int nump,double *fcp,double *dih_equ,double **fUMB);
double MD_Propagetor_NH_MP1998_AAFF_Amber_woH_UMB_wflag(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f,int bodflag,int angflag,int dihflag,int LJ14flag,int es14flag,int LJflag,int esflag,int *pairp,int nump,double *fcp,double *dih_equ,double **fUMB);

double MD_Propagetor_NH_MP1998_AAFF_Amber_woH_UMB_wflag(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f,int bodflag,int angflag,int dihflag,int LJ14flag,int es14flag,int LJflag,int esflag,int *pairp,int nump,double *fcp,double *dih_equ,double **fUMB);

double MD_Propagetor_NH_MP1998_AAFF_AMBER(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e,struct force *f);

double MD_Propagetor_NH_MP1998_AAFF_Amber_UMB(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f,int *pairp,int nump,double *fcp,double *dih_equ,double *pUMB,double **fUMB);

double MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008_b(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLMAA_PROTEINS2008 *e_GOLM);

double MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008_b_2013_08_31(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLMAA_PROTEINS2008 *e_GOLM,double fact_b,double fact_a,double fact_t,double fact_NC,double fact_NNC);

double MD_Propagetor_NH_MP1998_MB_FASYS(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3],double de, double d2, struct potential_GOLMAA_PROTEINS2008 *e_GOLM,int *dih_equ1,int *dih_equ2,double *e1,double *e2,double *p_MB,double **f_MB);

double MD_Propagetor_NH_MP1998_FASYS(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3],double de, double d2, struct potential_GOLMAA_PROTEINS2008 *e_GOLM,int *dih_equ1,int *dih_equ2,double *p_MB,double **f_MB,double *KE);

double MD_Propagetor_NH_MP1998_FASYS_2(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3],double de, double d2, struct potential_GOLMAA_PROTEINS2008 *e_GOLM,int *dih_equ1,int *dih_equ2,double *p_MB,double **f_MB);

#endif
