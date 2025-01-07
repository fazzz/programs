
#ifndef INCLUDE_TA_MD
#define INCLUDE_TA_MD

#include <netcdf.h>
#include "netcdf_mineL.h"

#include "FFL.h"

double runTACCM_CGAA_MD_NHC_MP1998_Amber_PROTEINS2008(// AA ////////////////////////////////////////////
						      double *crdAA,double *velAA, 
						      double *zetaAA,double *V_zetaAA, double QAA,
						      struct potential e, struct force f, double TAA, double NfKTAA,
						      double *avePEAA, double *aveKEAA,double *aveTAA,
						      double *varPEAA, double *varKEAA,double *varTAA,
						      struct my_netcdf_out_id_MCD nc_id_MCDAA,  FILE *outputfileAA,
						      // CG ////////////////////////////////////////////
						      double *crdCG,double *velCG, 
						      double *zetaCG,double *V_zetaCG, double QCG,
						      struct potential_GOLMAA_PROTEINS2008 e_CG, 
						      double TCG, double NfKTCG,
						      double *avePECG, double *aveKECG,double *aveTCG,
						      double *varPECG, double *varKECG,double *varTCG,
						      struct my_netcdf_out_id_MCD nc_id_MCDCG,  FILE *outputfileCG,
						      // Z  /////////////////////////////////////////////
						      double *Z,double *velZ,double massZ,
						      double *zetaZ,double *V_zetaZ,
						      double QZ,double NfKTZ,double TZ,
						      int numZ,double KZAA, double KZCG,int **pairs,
						      double *avePEZ, double *aveKEZ,double *aveTZ,
						      double *varPEZ, double *varKEZ,double *varTZ, 
						      FILE *trjfileZ, FILE *trjfilThetaAA, FILE *trjfilThetaCG,
						      // CM  ///////////////////////////////////////////////
						      double *mass, int numatom, int numheavyatom,
						      int numstep, int interval,int *l,
						      double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
						      double UNITT, double k_B,double pi,
						      double *PEZAA, double *PEZCG, double *PEZ);

double runTACCM_CGAA_MD_NHC_MP1998_Amber_PROTEINS2008_Z_asCA(// AA ////////////////////////////////////////////
							     double *crdAA,double *velAA, 
							     double *zetaAA,double *V_zetaAA, double QAA,
							     struct potential e, struct force f, double TAA, double NfKTAA,
							     double *avePEAA, double *aveKEAA,double *aveTAA,
							     double *varPEAA, double *varKEAA,double *varTAA,
							     struct my_netcdf_out_id_MCD nc_id_MCDAA,  FILE *outputfileAA,
							     // CG ////////////////////////////////////////////
							     double *crdCG,double *velCG, 
							     double *zetaCG,double *V_zetaCG, double QCG,
							     struct potential_GOLMAA_PROTEINS2008 e_CG, 
							     double TCG, double NfKTCG,
							     double *avePECG, double *aveKECG,double *aveTCG,
							     double *varPECG, double *varKECG,double *varTCG,
							     struct my_netcdf_out_id_MCD nc_id_MCDCG,  FILE *outputfileCG,
							     // Z  /////////////////////////////////////////////
							     double *Z,double *velZ,double massZ,
							     double *zetaZ,double *V_zetaZ,
							     double QZ,double NfKTZ,double TZ,
							     int numZ,double KZAA, double KZCG,int *index,
							     double *avePEZ, double *aveKEZ,double *aveTZ,
							     double *varPEZ, double *varKEZ,double *varTZ, 
							     struct my_netcdf_out_id_MCD trjfileZ, 
							     struct my_netcdf_out_id_MCD trjfilThetaAA,
							     struct my_netcdf_out_id_MCD trjfilThetaCG,
							     // CM  ///////////////////////////////////////////////
							     double *mass, int numatom, int numheavyatom,
							     int numstep, int interval,int *l,
							     double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
							     double UNITT, double k_B,double pi,
							     double *PEZAA, double *PEZCG, double *PEZ);

double runTACCM_CGAA_MD_NHC_MP1998_2PROTEINS2008(// CG1 ////////////////////////////////////////////
						 double *crdCG1,double *velCG1, 
						 double *zetaCG1,double *V_zetaCG1, double QCG1,
						 struct potential_GOLMAA_PROTEINS2008 e_CG1, 
						 double TCG1, double NfKTCG1,
						 double *avePECG1, double *aveKECG1,double *aveTCG1,
						 double *varPECG1, double *varKECG1,double *varTCG1,
						 struct my_netcdf_out_id_MCD nc_id_MCDCG1,  FILE *outputfileCG1,
						 // CG2 ////////////////////////////////////////////
						 double *crdCG2,double *velCG2, 
						 double *zetaCG2,double *V_zetaCG2, double QCG2,
						 struct potential_GOLMAA_PROTEINS2008 e_CG2, 
						 double TCG2, double NfKTCG2,
						 double *avePECG2, double *aveKECG2,double *aveTCG2,
						 double *varPECG2, double *varKECG2,double *varTCG2,
						 struct my_netcdf_out_id_MCD nc_id_MCDCG2,  FILE *outputfileCG2,
						 // Z  /////////////////////////////////////////////
						 double *Z,double *velZ,double massZ,
						 double *zetaZ,double *V_zetaZ,
						 double QZ,double NfKTZ,double TZ,
						 int numZ,double KZCG1, double KZCG2,int **pairs,
						 double *avePEZ, double *aveKEZ,double *aveTZ,
						 double *varPEZ, double *varKEZ,double *varTZ, 
						 FILE *trjfileZ, FILE *trjfilThetaCG1, FILE *trjfilThetaCG2,
						 // CM  ///////////////////////////////////////////////
						 double *mass, int numatom, int numheavyatom,
						 int numstep, int interval,int *l,
						 double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
						 double UNITT, double k_B,double pi,
						 double *PEZCG1, double *PEZCG2, double *PEZ);

#endif
