#ifndef INCLUDE_MD
#define INCLUDE_MD

#include "FFL.h"
#include "GOLMAA_hybrid_set.h"
#include "GOLMAA_hybrid.h"
#include "GOLM_Clementi_set.h"
#include "GOLM_Clementi.h"
#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"
#include "GOLMAA_MB_PROTEINS2008.h"

#define ON 1
#define OFF 0

#define NVE 0
#define NVT 1

int MODE;
int MODEV;

double GearsConstant[6];
double Telar_Matrix[6][6];

double MD_Generate_inivelo(double *vel,double *mass,int numatom,double KbT);

double MD_Propagetor_Iso(double *crd,double *vel,double *mass,int numatom,double IsoCoff,double dt,double *KE,double *PE);

double MD_Propagetor_Iso_JCP2003(double *crd,double *vel,double *mass,int numatom,double K,double dt,struct potential *e, struct force *f);

double MD_Propagetor_Iso_JCP2003_GOLMAA_JCTC2011(double *crd,double *vel,double *mass,int numatom,double K,double dt,struct potential *e, struct force *f, struct potential_GOLMAA_hybrid *e_GOLM);

double MD_Propagetor_Iso_JCP2003_GOLM_Clementi(double *crd,double *vel,double *mass,int numatom,double K,double dt, struct potential_GOLM_Clementi *e_GOLM);

double MD_Propagetor_vV_NVE(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f);

double MD_Propagetor_vV_NVE_GOLMAA_JCTC2011(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f, struct potential_GOLMAA_hybrid *e_GOLM);

double MD_Propagetor_vV_NVE_14LJdab(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f);

double MD_Propagetor_vV_NVE_GOLMAA_PROTEINS2008(double *crd,double *vel,double *mass,int numatom,double dt,struct potential_GOLMAA_PROTEINS2008 *e_GOLM);

double MD_Propagetor_vV_NVE_GOLMAA_JCTC2011_wonat(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f, struct potential_GOLMAA_hybrid *e_GOLM);

double MD_Propagetor_vV_NVE_14LJdba(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f);

double MD_Propagetor_vV_NVE_GOLM_Clementi(double *crd,double *vel,double *mass,int numatom,double dt,struct potential_GOLM_Clementi *e_GOLM);

double MD_Propagetor_vV_NVE_wc(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f,int bflag,int aflag,int dflag,int eflag, int LJflag,int e14flag,int LJ14flag );

double MD_Propagetor_vV_NVE_FASYS(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f);

double MD_Propagetor_vV_NVE_wflag(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f,int bodflag,int angflag,int dihflag,int LJ14flag,int es14flag,int LJflag,int esflag);

double MD_Propagetor_vV_NVE_woH_wflag(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f,int bodflag,int angflag,int dihflag,int LJ14flag,int es14flag,int LJflag,int esflag);

double MD_Propagetor_vV_NVE_GOLMAA_MB_PROTEINS2008(double *crd,double *vel,double *mass,int numatom,double dt,double de, double d2,struct potential_GOLMAA_MB_PROTEINS2008 *e_GOLM);

double MD_Propagetor_vV_NVE_AAFF_Amber(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f);

double MD_Propagetor_vV_NVE_UMB(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f,int *pairp,int nump,double *fcp,double *dih_equ,double *pUMB,double **fUMB);

double MD_Propagetor_vV_NVE_GOLMAA_PROTEINS2008_b(double *crd,double *vel,double *mass,int numatom,double dt,struct potential_GOLMAA_PROTEINS2008 *e_GOLM);

#endif

