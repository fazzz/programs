
#ifndef INCLUDE_TA_MD
#define INCLUDE_TA_MD

#include <netcdf.h>
#include "netcdf_mineL.h"

#include "FFL.h"

double runTACCM_MD_NHC_MP1998_Amber_AAFF(double *crd,double *vel, double *mass, int numatom,
					 double *zeta,double *V_zeta, double Q,
					 struct potential e, struct force f,
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

/*************************************************************************************************************************/
/* double runTACCM_MD_NHC_MP1998_GOLMAA_PROTEINS2008(double *crd,double *vel, double *mass, int numatom,		 */
/* 						  double *zeta,double *V_zeta, double Q,				 */
/* 						  struct potential e, struct force f,  					 */
/* 						  struct potential_GOLMAA_PROTEINS2008 e_GOLM,				 */
/* 						  double T, double NfKT, int numstep, int interval,int *l,		 */
/* 						  double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,		 */
/* 						  double *avePE, double *aveKE,double *aveT,				 */
/* 						  double *varPE, double *varKE,double *varT, 				 */
/* 						  double UNITT, double k_B,						 */
/* 						  struct my_netcdf_out_id_MCD nc_id_MCD,  FILE *outputfile,		 */
/* 						  //////////////// TACCM ///////////////////////			 */
/* 						  double *Z,double *velZ,double massZ,					 */
/* 						  double *zetaZ,double *V_zetaZ,					 */
/* 						  double TZ,double QZ,double NfKTZ,int numZ,				 */
/* 						  double Kapa,int **pairs,double pi,					 */
/* 						  double *avePEZ, double *aveKEZ,double *aveTZ,				 */
/* 						  double *varPEZ, double *varKEZ,double *varTZ, 			 */
/* 						  FILE *trjfileZ, FILE *trjfilTheta					 */
/* 						  //////////////// TACCM ///////////////////////			 */
/* 						  );									 */
/*************************************************************************************************************************/

#endif
