
#ifndef INCLUDE_TA_MD
#define INCLUDE_TA_MD

#include <netcdf.h>
#include "netcdf_mineL.h"

#include "FFL.h"

double runTACCM_CGAA_MD_NHC_MP1998_Amber_PROTEINS2008_Amber_hybrid(// AA ////////////////////////////////////////////
								   double *crdAA,double *velAA, 
								   double *zetaAA,double *V_zetaAA, double QAA,
								   struct potential e, struct force f,
								   double TAA, double NfKTAA,
								   double *avePEAA, double *aveKEAA,double *aveTAA,
								   double *varPEAA, double *varKEAA,double *varTAA,
								   struct my_netcdf_out_id_MCD nc_id_MCDAA,
								   FILE *outputfileAA,
								   // CG ////////////////////////////////////////////
								   double *crdCG,double *velCG, 
								   double *zetaCG,double *V_zetaCG, double QCG,
								   struct potential_GOLMAA_PROTEINS2008 e_CG, 
								   double TCG, double NfKTCG,
								   double *avePECG, double *aveKECG,double *aveTCG,
								   double *varPECG, double *varKECG,double *varTCG,
								   struct my_netcdf_out_id_MCD nc_id_MCDCG,
								   FILE *outputfileCG,
								   // Z  ///////////////////////////////////////////
								   double *Z,double *velZ,double massZ,
								   double *zetaZ,double *V_zetaZ,
								   double QZ,double NfKTZ,double TZ,
								   int numZ,double KZAA, double KZCG,int **pairs,
								   double *avePEZ, double *aveKEZ,double *aveTZ,
								   double *varPEZ, double *varKEZ,double *varTZ, 
								   FILE *trjfileZ, FILE *trjfilThetaAA,
								   FILE *trjfilThetaCG,
								   // CM  //////////////////////////////////////////
								   double *mass, int numatom, int numheavyatom,
								   int numstep, int interval,int *l,
								   double dt,double dt2,
								   double wdt2[3] ,double wdt4[3] ,int nc,
								   double UNITT, double k_B,double pi,
								   double *PEZAA, double *PEZCG, double *PEZ);

#endif
