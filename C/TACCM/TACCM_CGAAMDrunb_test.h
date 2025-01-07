
#ifndef INCLUDE_TA_MD
#define INCLUDE_TA_MD

#include <netcdf.h>
#include "netcdf_mineL.h"

#include "FFLc.h"

double runTACCM_CGAA_MDb_NHC_MP1998_test(// AA /////////////////////////////////////////////////////////
					 double *crdAA,double *velAA, 
					 double *zetaAA,double *V_zetaAA, double QAA,
					 struct potential e, struct force f, struct AmberParmL ap_AA,
					 double TAA, double NfKTAA,
					 double *avePEAA, double *aveKEAA,double *aveTAA,
					 double *varPEAA, double *varKEAA,double *varTAA,
					 struct my_netcdf_out_id_MCD nc_id_MCDAA,  FILE *outputfileAA,
					 // CG /////////////////////////////////////////////////////////
					 double *crdCG,double *velCG, 
					 double *zetaCG,double *V_zetaCG, double QCG,
					 struct potential e_CG, struct force f_CG, struct AmberParmL ap_CG,
					 double TCG, double NfKTCG,
					 double *avePECG, double *aveKECG,double *aveTCG,
					 double *varPECG, double *varKECG,double *varTCG,
					 struct my_netcdf_out_id_MCD nc_id_MCDCG,  FILE *outputfileCG,
					 // Z  /////////////////////////////////////////////////////////
					 double *Z,double *velZ,double massZ,
					 double *zetaZ,double *V_zetaZ,
					 double QZ,double NfKTZ,double TZ,
					 int numZ,double KZAA, double KZCG,int **pairs,
					 double *avePEZ, double *aveKEZ,double *aveTZ,
					 double *varPEZ, double *varKEZ,double *varTZ, 
					 FILE *trjfileZ, FILE *trjfilThetaAA, FILE *trjfilThetaCG,
					 // CM  /////////////////////////////////////////////////////////
					 double *mass, int numatom, int numstep, int interval,int *l,
					 double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
					 double UNITT, double k_B,double pi,
					 double *PEZAA, double *PEZCG, double *PEZ);

double runMDb_NHC_MP1998_Amber_AAFF_test(double *crd,double *vel, double *mass, int numatom,
					 double *zeta,double *V_zeta, double Q,
					 struct potential e, struct force f, struct AmberParmL ap,
					 double T, double NfKT, int numstep, int interval,int *l,
					 double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
					 double *avePE, double *aveKE,double *aveT,
					 double *varPE, double *varKE,double *varT, double UNITT, double k_B,
					 struct my_netcdf_out_id_MCD nc_id_MCD,  FILE *outputfile);

double runTACCM_MDb_NHC_MP1998_Amber_AAFF_test(double *crd,double *vel, double *mass, int numatom,
					       double *zeta,double *V_zeta, double Q,
					       struct potential e, struct force f, struct AmberParmL ap,
					       double T, double NfKT, int numstep, int interval,int *l,
					       double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
					       double *avePE, double *aveKE,double *aveT,
					       double *varPE, double *varKE,double *varT, double UNITT, double k_B,
					       struct my_netcdf_out_id_MCD nc_id_MCD,  FILE *outputfile,
					       //////////////// TACCM ///////////////////////
					       double *Z,double *velZ,double massZ,
					       double *zetaZ,double *V_zetaZ,
					       double TZ,double QZ,double NfKTZ,int numZ,
					       double Kapa,int **pairs,double pi,
					       double *avePEZ, double *aveKEZ,double *aveTZ,
					       double *varPEZ, double *varKEZ,double *varTZ, 
					       FILE *trjfileZ, FILE *trjfilTheta
					       //////////////// TACCM ///////////////////////
					       );

double CE_TACCMb_CGAA_test(double *crdAA,double *crdCG,double *Z, int numatom,int numZ,
			  double KZAA,double KZCG,int **pairs,double pi,
			   double *EAA,double *ECG,double *EZ);

double MDb_Propagetor_NH_MP1998_AAFF_Amber_test(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f, struct AmberParmL ap);

double TACCM_MDb_Propagetor_NH_MP1998_AA_test(double *crd,double *vel,double *mass,
					      double *zeta,double *V_zeta,double Q,
					      double NfKT,int numatom,double *KE,double *KEv,double *PEv,
					      double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
					      struct potential *e, struct force *f, struct AmberParmL ap,
					      double *Z,  int numZ,double *theta,double Kapa,
					      int **pairs, double *PEZ,double pi);

#endif
