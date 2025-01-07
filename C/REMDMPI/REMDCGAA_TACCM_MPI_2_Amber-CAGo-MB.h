

#ifndef INCLUDE_REMD_TAM
#define INCLUDE_REMD_TAM

#include <netcdf.h>
#include "mpi.h"

//#include "REMDMPI.h"
#include "REMD_TACCM_MPI.h"

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

#include "FFLc.h"
#include "PTLb.h"
#include "EF.h"

//#include "TOPO.h"
//#include "LA.h"

#include "GOLM_Clementi_MB.h"
#include "GOLM_Clementi_MB_wmCutOff.h"
#include "GOLM_Clementi_set_wmCutOff.h"

#include "RAND.h"
#include "BOXMULL.h"
#include "MC.h"

#include "MDrunb.h"
#include "TACCM.h"
#include "TACCM_CGAAMDrun_Amber_CAGo_MB.h"
//#include "TACCM_CGAAMDrun_Amber_CAGo_MB_2.h"
#include "TACCM_CGAAMDrun_Amber_CAGo_MB_CTheta.h"

#include "netcdf_mineL.h"

#define AAINPF 0
#define CGINPF 1
#define AAKZ 2
#define CGKZ 3

struct AADataforREMD_Amber_CAGo_MB{  
  double *crd,*vel;
  double zeta,V_zeta,Q;
  struct potential e;
  struct force f;
  double T,NfKT;
  double *avePE,*aveKE,*aveT;
  double *varPE,*varKE,*varT;

  struct my_netcdf_out_id_MCD nc_id_MCD;  
  FILE *outputfile;
};

struct CGDataforREMD_Amber_CAGo_MB{  
  double *crd,*vel;
  double zeta,V_zeta,Q;
  struct potential e;
  struct force f;
  struct potential_GOLM_Clementi_MB e_GOLM;
  double T,NfKT;
  double *avePE,*aveKE,*aveT;
  double *varPE,*varKE,*varT;

  struct my_netcdf_out_id_MCD nc_id_MCD;  
  FILE *outputfile;
};

struct AACGCommonDataforREMD_Amber_CAGo_MB{  
  int numatom,numheavyatom,numCAatom,numres;

  double *mass,*massCA;  
  int numstep,interval;
};

struct TACCMDataforREMD_Amber_CAGo_MB{  
  int numZ;

  int numZ_dih,numZ_ang,numZ_bon;

  double *Z,*velZ,massZ;
  double zetaZ,V_zetaZ;

  double T,QZ,NfKTZ;
  double KZAA,KZCG;  
  int **pairs;
  int **pairs_dihed_AA,**pairs_angle_AA,**pairs_bond_AA;
  int **pairs_dihed_CG,**pairs_angle_CG,**pairs_bond_CG;

  double *avePEZ,*aveKEZ,*aveTZ;
  double *varPEZ,*varKEZ,*varTZ;

  FILE *trjfileZ,*trjfilThetaAA,*trjfilThetaCG;
};

double **MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_Amber_CAGo_MB(int myrank,int num_procs,int tag, MPI_Status* status,
							    int numRE, int numEX, double *KZAA, double *KZCG,
							    struct AADataforREMD_Amber_CAGo_MB AAData,
							    struct CGDataforREMD_Amber_CAGo_MB CGData,
							    struct TACCMDataforREMD_Amber_CAGo_MB ZData,
							    struct AACGCommonDataforREMD_Amber_CAGo_MB CData,
							    struct AmberParmL ap_AA, struct AmberParmL ap_CG,
							    double de, double d2,
							    double T0AA,double T0CG, double T0Z, 
							    int numstep, int interval, 
							    double dt,double dt2,
							    double wdt2[3],double wdt4[3], int nc,
							    double UNITT, double k_B, double tau, double pi,
							    double parameterCG, FILE* logfile );


void CGAAREMDreadInputs_Amber_CAGo_MB(FILE *inputfile,int numatom,int numRE,int myrank,
				      double *crdAA,double *velAA, double *crdCG,double *velCG,
				      double *KZAA, double *KZCG, struct AmberParmL ap_CG);

/****************************************************************************************************************************/
/* double runTACCM_CGAA_MD_NHC_MP1998_Amber_CAGo_MB(// AA /////////////////////////////////////////////////////////	    */
/* 						 double *crdAA,double *velAA, 						    */
/* 						 double *zetaAA,double *V_zetaAA, double QAA,				    */
/* 						 struct potential e, struct force f, double TAA, double NfKTAA,		    */
/* 						 double *avePEAA, double *aveKEAA,double *aveTAA,			    */
/* 						 double *varPEAA, double *varKEAA,double *varTAA,			    */
/* 						 struct my_netcdf_out_id_MCD nc_id_MCDAA,  FILE *outputfileAA,		    */
/* 						 // CG /////////////////////////////////////////////////////////	    */
/* 						 double *crdCG,double *velCG, 						    */
/* 						 double *zetaCG,double *V_zetaCG, double QCG,				    */
/* 						 struct potential_GOLM_Clementi_MB e_CG,				    */
/* 						 double de, double d2,							    */
/* 						 double TCG, double NfKTCG,						    */
/* 						 double *avePECG, double *aveKECG,double *aveTCG,			    */
/* 						 double *varPECG, double *varKECG,double *varTCG,			    */
/* 						 struct my_netcdf_out_id_MCD nc_id_MCDCG,  FILE *outputfileCG,		    */
/* 						 // Z  /////////////////////////////////////////////////////////	    */
/* 						 double *Z,double *velZ,double massZ,					    */
/* 						 double *zetaZ,double *V_zetaZ,						    */
/* 						 double QZ,double NfKTZ,double TZ,					    */
/* 						 double KZAA, double KZCG,						    */
/* 						 int numZ_dih,int **pairs_dihe_AA,int **pairs_dihe_CG,			    */
/* 						 int numZ_ang,int **pairs_angl_AA,int **pairs_angl_CG,			    */
/* 						 int numZ_bon,int **pairs_bond_AA,int **pairs_bond_CG,			    */
/* 						 double *avePEZ, double *aveKEZ,double *aveTZ,				    */
/* 						 double *varPEZ, double *varKEZ,double *varTZ, 				    */
/* 						 FILE *trjfileZ, FILE *trjfilThetaAA, FILE *trjfilThetaCG,		    */
/* 						 // CM  /////////////////////////////////////////////////////////	    */
/* 						 double *mass,double *massCA,						    */
/* 						 int numatom, int numCAatom, 						    */
/* 						 int numstep, int interval,int *l,					    */
/* 						 double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,		    */
/* 						 double UNITT, double k_B,double pi,					    */
/* 						 double *PEZAA, double *PEZCG, double *PEZ);				    */
/****************************************************************************************************************************/

/*******************************************************************************************************************************/
/* double TACCM_MD_Propagetor_NH_MP1998_AA_Amber_CAGo_MB(double *crd,double *vel,double *mass,				       */
/* 						      double *zeta,double *V_zeta,double Q,				       */
/* 						      double NfKT,int numatom,double *KE,double *KEv,double *PEv,	       */
/* 						      double dt,double dt2,int nc,double wdt4[3],double wdt2[3],	       */
/* 						      struct potential *e, struct force *f,				       */
/* 						      double *Z, int numZ_dih,int numZ_ang, int numZ_bon,		       */
/* 						      double *theta,double Kapa,					       */
/* 						      int **pairs_dihe_AA, int **pairs_angl_AA, int **pairs_bond_AA,	       */
/* 						      double *PEZ,double pi);						       */
/* 															       */
/* double TACCM_MD_Propagetor_NH_MP1998_CG_Amber_CAGo_MB(double *crd,double *vel,double *mass,				       */
/* 						      double *zeta,double *V_zeta,double Q,				       */
/* 						      double NfKT,int numCAatom,double *KE,double *KEv,double *PEv,	       */
/* 						      double dt,double dt2,int nc,double wdt4[3],double wdt2[3],	       */
/* 						      struct potential_GOLM_Clementi_MB *e_CG,				       */
/* 						      double de, double d2,						       */
/* 						      double *Z, int numZ_dih,int numZ_ang, int numZ_bon,		       */
/* 						      double *theta,double Kapa,					       */
/* 						      int **pairs_dihe_CG, int **pairs_angl_CG, int **pairs_bond_CG,	       */
/* 						      double *PEZ,double pi);						       */
/* 															       */
/* double TACCM_CTheta_Amber_CAGo_MB(double *crd,int numatom,double *theta, 						       */
/* 				  int numdihe, int **pairs_dih_AA,							       */
/* 				  int numangl, int **pairs_ang_AA,							       */
/* 				  int numbond, int **pairs_bon_AA, 							       */
/* 				  double pi);										       */
/* 															       */
/* double TACCM_calc_eff_FF_MD_Amber_CAGo_MB(double *crd,int numatom, double *theta, double *Z,  			       */
/* 					  int numZ_dih,int numZ_ang,int numZ_bon,					       */
/* 					  double Kapa, double **frcZ,							       */
/* 					  int **pairs_dih,int **pairs_ang,int **pairs_bon,				       */
/* 					  double pi);									       */
/*******************************************************************************************************************************/

/****************************************************************************************************************/
/* double CE_TACCMb_CGAA_Amber_CAGo_MB(double *crdAA,double *crdCG,double *Z, int numatom,int numCAatom,        */
/* 				    int numZ_dih, int **pairs_dihed_AA, int **pairs_dihed_CG,		        */
/* 				    int numZ_ang, int **pairs_angle_AA, int **pairs_angle_CG,		        */
/* 				    int numZ_bon, int **pairs_bond_AA,  int **pairs_bond_CG,		        */
/* 				    double KZAA,double KZCG, double pi,					        */
/* 				    double *EAA,double *ECG,double *EZ);				        */
/****************************************************************************************************************/

#endif

