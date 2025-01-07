

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
#include "TACCM_CGAAMDrun_Amber_PROTEINS2008_Amber_hybrid.h"
//#include "TACCM_CGAAMDrun_test_CG.h"

#include "netcdf_mineL.h"


double **MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_Amber_PROTEINS2008_Amber_hybrid(int myrank,int num_procs,
									      int tag, MPI_Status* status,
									      int numRE, int numEX, 
									      double *KZAA, double *KZCG,
									      struct AADataforREMD_Amber AAData,
									      struct CGDataforREMD_PROTEINS2008 CGData,
									      struct TACCMDataforREMD_A_P2008 ZData,
									      struct AACGCommonDataforREMD_A_P2008 CData,
									      double T0AA,double T0CG, double T0Z, 
									      int numstep, int interval, 
									      double dt,double dt2,
									      double wdt2[3],double wdt4[3], int nc,
									      double UNITT, double k_B, double tau,
									      double pi, FILE* logfile );

void  CGAAREMDreadInputs_Amber_PROTEINS2008_Amber_hybrid(FILE *inputfile,int numatom,int numRE,int myrank,
							 double *crdAA,double *velAA, double *crdCG,double *velCG,
							 double *KZAA, double *KZCG);

#endif
