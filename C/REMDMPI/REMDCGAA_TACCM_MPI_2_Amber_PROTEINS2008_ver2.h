

#ifndef INCLUDE_REMD_TAM
#define INCLUDE_REMD_TAM

#include <netcdf.h>
#include "mpi.h"

#include "REMDMPI.h"
#include "REMD_TACCM_MPI.h"

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

#include "FFL.h"
#include "PTL.h"
#include "EF.h"

#include "RAND.h"
#include "BOXMULL.h"
#include "MC.h"

#include "MDrun.h"
#include "TACCM.h"
#include "TACCM_CGAAMDrun_Amber_PROTEINS2008.h"
//#include "TACCM_CGAAMDrun_test_CG.h"

#include "netcdf_mineL.h"

#include "ABAb.h"

/***********************************************************************************************/
/* struct TACCMDataforREMD_A_P2008_ver2{  						       */
/*   int numZ;										       */
/*   double *Z,massZ;									       */
/* 											       */
/*   double KZAA,KZCG;									       */
/*   int **pairs;									       */
/* 											       */
/*   double *crd,*vel;									       */
/*   double zeta,V_zeta,Q;								       */
/*   struct potential e;								       */
/*   struct force f;									       */
/*   double T,NfKT;									       */
/*   double *avePE,*aveKE,*aveT;							       */
/*   double *varPE,*varKE,*varT;							       */
/* 											       */
/*   CLTb *clt;										       */
/*   double *q,*qvel;									       */
/*   double *predict,*correct;								       */
/*   double s,s_vel,s_acc,gzi,gzi_vel,predict_gzi[5],correct_gzi[5],predict_s[5],correct_s[5]; */
/*   double *vel_Term,**predict_Term,**predict_Term2,**correct_Term,**correct_Term2;	       */
/*   double To,KEo;									       */
/* 											       */
/*   struct my_netcdf_out_id_MCD nc_id_MCD;  FILE *outputfile;				       */
/* 											       */
/* };											       */
/***********************************************************************************************/

double **MPI_CGFGTREM_TACCM_ABAbMD_NH_new_Amber_PROTEINS2008_ver2(int myrank,int num_procs,
								  int tag, MPI_Status* status,
								  int numRE, int numEX, double *KZAA, double *KZCG,
								  struct AADataforREMD_Amber AAData,
								  struct CGDataforREMD_PROTEINS2008 CGData,
								  struct TACCMDataforREMD_A_P2008_ver2 ZData,
								  struct AACGCommonDataforREMD_A_P2008 CData,
								  double T0AA,double T0CG, double T0Z, 
								  int numstep, int interval, 
								  double dt,double tau, double tau2,
								  double UNITT, double k_B, 
								  double pi, FILE* logfile );

#endif
