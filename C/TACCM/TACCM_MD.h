

#ifndef INCLUDE_FVDIHED
#define INCLUDE_FVDIHED

double TACCM_calc_eff_FF_MD(double *crd,int numatom, double *theta, double *Z,  int numZ,double Kapa, double **frcZ, int **pairs,double pi);

double TACCM_MD_Propagetor_NH_MP1998_AAFF_Amber(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KE,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f,double *Z,  int numZ,double *theta,double Kapa,int **pairs, double *PEZ,double pi);

double TACCM_MD_Propagetor_NH_MP1998_Z(double *Z,double *velZ,double massZ,double *theta,
				       double *zeta,double *V_zeta,double Q,double NfKT,int numZ,
				       double *KE, double *KEv,double *PEv,
				       double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
				       double KZ,double *PEZ,double *f,double pi);

double TACCM_MD_Propagetor_vV_NVE_Z(double *Z, double *velZ,double massZ,double *theta,int numZ,
				    double dt, double KZ,double *PEZ,double *f,double pi);

double TACCM_MD_Propagetor_NH_Single_part_MP1996(double *vel,double massZ,double *zeta,double *V_zeta,double Q,double NfKT,int numZ,int nc,double wdt4[3],double wdt2[3]);


double TACCM_MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KE,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLMAA_PROTEINS2008 *e_GOLM,double *Z,  int numZ,double *theta,double Kapa,int **pairs, double *PEZ,double pi);

double TACCM_MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008_2013_08_31(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KE,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLMAA_PROTEINS2008 *e_GOLM,double *Z,  int numZ,double *theta,double Kapa,int **pairs, double *PEZ,double pi, double fact_b,double fact_a,double fact_t,double fact_NC,double fact_NNC);

double TACCM_MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008_Amber_hybrid(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KE,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLMAA_PROTEINS2008 *e_GOLM,double *Z,  int numZ,double *theta,double Kapa,int **pairs, double *PEZ,double pi);

double TACCM_MD_Propagetor_NH_MP1998_AAFF_Amber_Z_asCA(double *crd,double *vel,double *mass,
						       double *zeta,double *V_zeta,
						       double Q,double NfKT,int numatom,
						       double *KE,double *KEv,double *PEv,
						       double dt,double dt2,int nc,double wdt4[3],double wdt2[3], 
						       struct potential *e, struct force *f,
						       double *Z,  int numZ,double *theta,
						       double Kapa,int *index, double *PEZ,double pi);

double TACCM_MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008_Z_asCA(double *crd,double *vel,double *mass,
								double *zeta,double *V_zeta,double Q,double NfKT,
								int numatom,double *KE,double *KEv,double *PEv,
								double dt,double dt2,int nc,
								double wdt4[3],double wdt2[3], 
								struct potential_GOLMAA_PROTEINS2008 *e_GOLM,
								double *Z,  int numZ,double *theta,
								double Kapa,int *index, double *PEZ,double pi);

#endif
