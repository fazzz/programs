
#ifndef INCLUDE_TA_MD
#define INCLUDE_TA_MD

#include <netcdf.h>
#include "netcdf_mineL.h"

#include "FFL.h"
#include "GOLM_Clementi_MB.h"
#include "GOLM_Clementi_MB_wmCutOff.h"
#include "GOLM_Clementi_set_wmCutOff.h"

double runTACCM_CGAA_MD_NHC_MP1998_Amber_CAGo_MB(// AA /////////////////////////////////////////////////////////
						 double *crdAA,double *velAA, 
						 double *zetaAA,double *V_zetaAA, double QAA,
						 struct potential e, struct force f, double TAA, double NfKTAA,
						 double *avePEAA, double *aveKEAA,double *aveTAA,
						 double *varPEAA, double *varKEAA,double *varTAA,
						 struct my_netcdf_out_id_MCD nc_id_MCDAA,  FILE *outputfileAA,
						 // CG /////////////////////////////////////////////////////////
						 double *crdCG,double *velCG, 
						 double *zetaCG,double *V_zetaCG, double QCG,
						 struct potential_GOLM_Clementi_MB e_CG,
						 double de, double d2,
						 double TCG, double NfKTCG,
						 double *avePECG, double *aveKECG,double *aveTCG,
						 double *varPECG, double *varKECG,double *varTCG,
						 struct my_netcdf_out_id_MCD nc_id_MCDCG,  FILE *outputfileCG,
						 // Z  /////////////////////////////////////////////////////////
						 double *Z,double *velZ,double massZ,
						 double *zetaZ,double *V_zetaZ,
						 double QZ,double NfKTZ,double TZ,
						 double KZAA, double KZCG,
						 int numZ_dih,int **pairs_dihe_AA,int **pairs_dihe_CG,
						 int numZ_ang,int **pairs_angl_AA,int **pairs_angl_CG,
						 int numZ_bon,int **pairs_bond_AA,int **pairs_bond_CG,
						 double *avePEZ, double *aveKEZ,double *aveTZ,
						 double *varPEZ, double *varKEZ,double *varTZ, 
						 FILE *trjfileZ, FILE *trjfilThetaAA, FILE *trjfilThetaCG,
						 // CM  /////////////////////////////////////////////////////////
						 struct AmberParmL ap_AA,struct AmberParmL ap_CG,
						 double *mass,double *massCA,
						 int numatom, int numCAatom, 
						 int numstep, int interval,int *l,
						 double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
						 double UNITT, double k_B,double pi,
						 double *PEZAA, double *PEZCG, double *PEZ);

double TACCM_MD_Propagetor_NH_MP1998_AA_Amber_CAGo_MB(double *crd,double *vel,double *mass,
						      double *zeta,double *V_zeta,double Q,
						      double NfKT,int numatom,double *KE,double *KEv,double *PEv,
						      double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
						      struct potential *e, struct force *f,
						      double *Z, int numZ_dih,int numZ_ang, int numZ_bon,
						      double *theta,double Kapa,
						      int **pairs_dihe_AA, int **pairs_angl_AA, int **pairs_bond_AA,
						      double *PEZ,double pi);

double TACCM_MD_Propagetor_NH_MP1998_CG_Amber_CAGo_MB(double *crd,double *vel,double *mass,
						      double *zeta,double *V_zeta,double Q,
						      double NfKT,int numCAatom,double *KE,double *KEv,double *PEv,
						      double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
						      struct potential_GOLM_Clementi_MB *e_CG,
						      double de, double d2,
						      double *Z, int numZ_dih,int numZ_ang, int numZ_bon,
						      double *theta,double Kapa,
						      int **pairs_dihe_CG, int **pairs_angl_CG, int **pairs_bond_CG,
						      double *PEZ,double pi);

double TACCM_calc_eff_FF_MD_Amber_CAGo_MB(double *crd,int numatom, double *theta, double *Z,  
					  int numZ_dih,int numZ_ang,int numZ_bon,
					  double Kapa, double **frcZ,
					  int **pairs_dih,int **pairs_ang,int **pairs_bon,
					  double pi);

double CE_TACCMb_CGAA_Amber_CAGo_MB(double *crdAA,double *crdCG,double *Z, int numatom,int numCAatom,
				    //int numZ_dih, int **pairs_dihed_AA, int **pairs_dihed_CG,
				    //int numZ_ang, int **pairs_angle_AA, int **pairs_angle_CG,
				    //int numZ_bon, int **pairs_bond_AA,  int **pairs_bond_CG,
				    int numZ_dih, int numZ_ang, int numZ_bon, 
				    int **pairs_dihed_AA, int **pairs_dihed_CG,
				    int **pairs_angle_AA, int **pairs_angle_CG,
				    int **pairs_bond_AA,  int **pairs_bond_CG,
				    double KZAA,double KZCG, double pi,
				    double *EAA,double *ECG,double *EZ);

/************************************************************************************/
/* double TACCM_CTheta_Amber_CAGo_MB(double *crd,int numatom,double *theta, 	    */
/* 				  int numdihe, int **pairs_dih_AA,		    */
/* 				  int numangl, int **pairs_ang_AA,		    */
/* 				  int numbond, int **pairs_bon_AA, 		    */
/* 				  double pi);					    */
/************************************************************************************/

double runTACCM_CGAA_MD_NHC_MP1998_Amber_CAGo_MB_2(// AA /////////////////////////////////////////////////////////
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
						   struct potential_GOLM_Clementi_MB e_CG,
						   double de, double d2,
						   struct AmberParmL ap_CG,
						   double TCG, double NfKTCG,
						   double *avePECG, double *aveKECG,double *aveTCG,
						   double *varPECG, double *varKECG,double *varTCG,
						   struct my_netcdf_out_id_MCD nc_id_MCDCG,  FILE *outputfileCG,
						   // Z  /////////////////////////////////////////////////////////
						   double *Z,double *velZ,double massZ,
						   double *zetaZ,double *V_zetaZ,
						   double QZ,double NfKTZ,double TZ,
						   int numZ_dih,int numZ_ang,int numZ_bon,
						   double KZAA, double KZCG,
						   int **pairs_dihed_AA,int **pairs_dihed_CG,
						   int **pairs_angle_AA,int **pairs_angle_CG,
						   int **pairs_bond_AA, int **pairs_bond_CG,
						   double *avePEZ, double *aveKEZ,double *aveTZ,
						   double *varPEZ, double *varKEZ,double *varTZ, 
						   FILE *trjfileZ, FILE *trjfilThetaAA, FILE *trjfilThetaCG,
						   // CM  /////////////////////////////////////////////////////////
						   double *mass, double *massCA, int numatom, int numCAatom,
						   int numstep, int interval,int *l,
						   double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
						   double UNITT, double k_B,double pi,
						   double *PEZAA, double *PEZCG, double *PEZ,
						   double *ele_ele,    
						   double *ALJ_s,    double *BLJ_s,
						   double *ele_ele_14, 
						   double *ALJ_14_s, double *BLJ_14_s);

double TACCM_CTheta_Amber_CAGo_MB_2(double *crd,int numatom,double *theta, 
				    int numdihe, int **pairs_dih_AA,
				    int numangl, int **pairs_ang_AA,
				    int numbond, int **pairs_bon_AA, 
				    double pi);

double TACCM_MD_Propagetor_NH_MP1998_AA_Amber_CAGo_MB_2(double *crd,double *vel,double *mass,
							double *zeta,double *V_zeta,double Q,
							double NfKT,int numatom,double *KE,double *KEv,double *PEv,
							double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
							struct potential *e, struct force *f, struct AmberParmL ap,
							double *Z,  int numZ_dih, int numZ_ang, int numZ_bon,
							double *theta,double Kapa,
							int **pairs_dih, int **pairs_ang, int **pairs_bon,
							double *PEZ,double pi,
							double *ele_ele,    double *ALJ_s,    double *BLJ_s,
							double *ele_ele_14, double *ALJ_14_s, double *BLJ_14_s);

double TACCM_MD_Propagetor_NH_MP1998_CG_Amber_CAGo_MB_2(double *crd,double *vel,double *mass,
							double *zeta,double *V_zeta,double Q,
							double NfKT,int numCAatom,
							double *KE,double *KEv,double *PEv,
							double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
							struct potential_GOLM_Clementi_MB *e_CG,
							double de, double d2,
							double *Z, int numZ_dih, int numZ_ang, int numZ_bon,
							double *theta,double Kapa,
							int **pairs_dih, int **pairs_ang, int **pairs_bon, 
							double *PEZ,double pi);

double TACCM_CGAA_MD_Propagetor_NH_MP1998_Z_Amber_CAGo_MB_2(double *Z,double *velZ,double massZ,
							    double *thetaAA,double *thetaCG,
							    double *zeta,double *V_zeta,double Q,double NfKT,
							    int numZ_dih, int numZ_ang, int numZ_bon,
							    double *KE, double *KEv,double *PEv,
							    double dt,double dt2,int nc,
							    double wdt4[3],double wdt2[3],
							    double KZAA, double KZCG,
							    double *PEZAA,double *PEZCG,double *PEZ,double *f,
							    double UNITT, double pi);

double TACCM_calc_eff_FF_MD_Amber_CAGo_MB(double *crd,int numatom, double *theta, double *Z,  
					  int numZ_dih,int numZ_ang,int numZ_bon,
					  double Kapa, double **frcZ,
					  int **pairs_dih,int **pairs_ang,int **pairs_bon,
					  double pi);

double TACCM_calc_eff_FF_Z_2_Amber_CAGo_MB(double *Z,int numZ_dih, int numZ_ang, int numZ_bon,
					   double *theta,double KZ,double *f, double *PE,
					   double UNITT, double pi);

#endif
