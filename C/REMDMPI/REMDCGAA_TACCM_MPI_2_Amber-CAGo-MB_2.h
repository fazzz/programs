

#ifndef INCLUDE_REMD_TAM
#define INCLUDE_REMD_TAM

#include <netcdf.h>
#include "mpi.h"

#include "REMDMPI.h"
#include "REMD_TACCM_MPI.h"

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

#include "FFLc.h"
#include "PTLb.h"
#include "EF.h"

#include "RAND.h"
#include "BOXMULL.h"
#include "MC.h"

#include "MDrunb.h"
#include "TACCM.h"
#include "TACCM_CGAAMDrunb_test.h"
//#include "TACCM_CGAAMDrun_Amber_CAGo_MB_2.h"

#include "netcdf_mineL.h"

#include "REMDCGAA_TACCM_MPI_2_testb.h"
#include "REMDCGAA_TACCM_MPI_2_Amber-CAGo-MB.h"

double CE_TACCMb_CGAA_Amber_CAGo_MB_2(double *crdAA,double *crdCG,double *Z,
				      int numatom,int numCAatom, int numZ,
				      double KZAA,double KZCG,
				      int **pairs_AA, int **pairs_CG,
				      double pi, double *EAA,double *ECG,double *EZ);

double CE_TACCMb_CGAA_Amber_CAGo_MB_3(double *crdAA,double *crdCG,double *Z,
				      int numatom,int numCAatom, 
				      int numZ_dih, int numZ_ang, int numZ_bon,
				      double KZAA,double KZCG,
				      int **pairs_dihed_AA, int **pairs_dihed_CG,
				      int **pairs_angle_AA, int **pairs_angle_CG,
				      int **pairs_bond_AA,  int **pairs_bond_CG,
				      double pi, double *EAA,double *ECG,double *EZ){

double **MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_Amber_CAGo_MB_2(int myrank,int num_procs,int tag, MPI_Status* status,
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
							      double parameterCG, FILE* logfile,
							      double *ele_ele,    
							      double *ALJ_s,    double *BLJ_s,
							      double *ele_ele_14, 
							      double *ALJ_14_s, double *BLJ_14_s);

#endif
