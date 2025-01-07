

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

#include "netcdf_mineL.h"

#include "REMDCGAA_TACCM_MPI_2_testb.h"

double **MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_testc(int myrank,int num_procs,int tag, MPI_Status* status,
						    int numRE, int numEX, double *KZAA, double *KZCG,
						    struct AADataforREMD_testb AAData,
						    struct AADataforREMD_testb CGData,
						    struct TACCMDataforREMD_testb ZData,
						    struct AACGCommonDataforREMD_testb CData,
						    struct AmberParmL ap_AA, struct AmberParmL ap_CG,
						    double T0AA,double T0CG, double T0Z, int numstep, int interval, 
						    double dt,double dt2,double wdt2[3],double wdt4[3], int nc,
						    double UNITT, double k_B, double tau, double pi,
						    double parameterCG, FILE* logfile );

void  CGAAREMDreadInputs_testc(FILE *inputfile,int numatom,int numRE,int myrank,
			      double *crdAA,double *velAA, double *crdCG,double *velCG,
			      double *KZAA, double *KZCG);

#endif
