
#ifndef INCLUDE_TA_MD
#define INCLUDE_TA_MD

#include <netcdf.h>
#include "netcdf_mineL.h"

#include "FFL.h"

double runTACCM_CGAA_MD_NHC_MP1998(// AA /////////////////////////////////////////////////////////
				   double *crdAA,double *velAA, 
				   double *zetaAA,double *V_zetaAA, double QAA,
				   struct potential e, struct force f, double TAA, double NfKTAA,
				   double *avePEAA, double *aveKEAA,double *aveTAA,
				   double *varPEAA, double *varKEAA,double *varTAA,
				   struct my_netcdf_out_id_MCD nc_id_MCDAA,  FILE *outputfileAA,
				   // CG /////////////////////////////////////////////////////////
				   double *crdCG,double *velCG, 
				   double *zetaCG,double *V_zetaCG, double QCG,
				   struct potential_GOLMAA_PROTEINS2008 e_GOLM, double TCG, double NfKTCG,
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

double TACCM_CGAA_MD_Propagetor_NH_MP1998_Z(double *Z,double *velZ,double massZ,double *thetaAA,double *thetaCG,
					    double *zeta,double *V_zeta,double Q,double NfKT,int numZ,
					    double *KE, double *KEv,double *PEv,
					    double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
					    double KZAA, double KZCG,double *PEZ,double *f,double pi);

double TACCM_CGAA_MD_Propagetor_NH_MP1998_Z_2(double *Z,double *velZ,double massZ,double *thetaAA,double *thetaCG,
					      double *zeta,double *V_zeta,double Q,double NfKT,int numZ,
					      double *KE, double *KEv,double *PEv,
					      double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
					      double KZAA, double KZCG,
					      double *PEAA,double *PECG,double *PEZ,double *f,double pi);

double CE_TACCM_CGAA(double *crdAA,double *crdCG,double *Z, int numatom,int numZ,
		     double KZAA,double KZCG,int **pairs,double pi,
		     double *EAA,double *ECG,double *EZ);

double TACCM_CGAA_MD_Propagetor_NH_MP1998_Z_2_asCA(double *Z,double *velZ,double massZ,
						   double *thetaAA,double *thetaCG,
						   double *zeta,double *V_zeta,double Q,double NfKT,int numZ,
						   double *KE, double *KEv,double *PEv,
						   double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
						   double KZAA, double KZCG,
						   double *PEZAA,double *PEZCG,double *PEZ,double *f,double pi);

#endif
