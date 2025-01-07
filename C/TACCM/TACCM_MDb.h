

#ifndef INCLUDE_FVDIHEDb
#define INCLUDE_FVDIHEDb

double TACCMb_calc_eff_FF_MD(double *crd,int numatom, double *theta, double *Z,  int numZ,double Kapa, double **frcZ, int **pairs,double pi);

double TACCMb_MD_Propagetor_NH_MP1998_AAFF_Amber(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KE,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f,double *Z,  int numZ,double *theta,double Kapa,int **pairs, double *PEZ,double pi, struct AmberParmL ap);

double TACCMb_MD_Propagetor_NH_MP1998_Z(double *Z,double *velZ,double massZ,double *theta,
				       double *zeta,double *V_zeta,double Q,double NfKT,int numZ,
				       double *KE, double *KEv,double *PEv,
				       double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
				       double KZ,double *PEZ,double *f,double pi);

double TACCMb_MD_Propagetor_vV_NVE_Z(double *Z, double *velZ,double massZ,double *theta,int numZ,
				    double dt, double KZ,double *PEZ,double *f,double pi);

double TACCMb_MD_Propagetor_NH_Single_part_MP1996(double *vel,double massZ,double *zeta,double *V_zeta,double Q,double NfKT,int numZ,int nc,double wdt4[3],double wdt2[3]);


double TACCMb_MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KE,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLMAA_PROTEINS2008 *e_GOLM,double *Z,  int numZ,double *theta,double Kapa,int **pairs, double *PEZ,double pi);

#endif
