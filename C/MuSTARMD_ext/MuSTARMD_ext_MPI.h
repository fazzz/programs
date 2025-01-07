
#ifndef INCLUDE_REMD_TAM
#define INCLUDE_REMD_TAM

#include <netcdf.h>
#include "mpi.h"

#include "REMDMPI.h"
#include "REMD_TACCM_MPI.h"

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

#define AAKZ 2
#define CGKZ 3

struct AAData_MuSTARMD{  
  double *crd,*vel;
  double zeta,V_zeta,Q;
  struct potential e;
  struct force f;
  double T,NfKT;
  double *avePE,*aveKE,*aveT;
  double *varPE,*varKE,*varT;

  struct my_netcdf_out_id_MCD nc_id_MCD;  FILE *outputfile;
};

struct AACGCommonData_MuSTARMD{  
  int numatom,numheavyatom,numres;

  double *mass;  int numstep,interval;
};

struct TACCMData_MuSTARMD{  
  int numZ;
  double *Z,*velZ,massZ;
  double zetaZ,V_zetaZ;

  double T,QZ,NfKTZ;
  double *KZAlpha;  int **pairs;

  double *avePEZ,*aveKEZ,*aveTZ;
  double *varPEZ,*varKEZ,*varTZ;

  FILE *trjfileZ,**trjfilThetaAlpha;
};

double ***MPI_MuSTARMD_ext_pep_NHC_MP1998(int ND,int *myrank,
					  int num_procs, int tag, MPI_Status* status,
					  int *numRE, int numEX,
					  double **KZAlpha,
					  struct AAData_MuSTARMD *AlphaData,
					  struct TACCMData_MuSTARMD ZData,
					  struct AACGCommonData_MuSTARMD CData,
					  struct AmberParmL *ap_Alpha,
					  double *T0Alpha, double T0Z, int numstep, int interval, 
					  double dt,double dt2,double wdt2[3],double wdt4[3], int nc,
					  double UNITT, double k_B, double tau, double pi,
					  FILE* logfile );

void  MuSTARMD_ext_readInputs(int ND, FILE *inputfile,
			      int numatom,int *numRE,int myrank,
			      double **crdAlpha,double **velAlpha, double **KZAlpha);

double runMuSTARMD_AmberType_NHC_MP1998(int ND,
					struct AAData_MuSTARMD *AlphaData,
					struct TACCMData_MuSTARMD ZData,
					struct AACGCommonData_MuSTARMD CData,
					int numstep,int interval,int l,
					double dt,double dt2,double wdt2[3],double wdt4[3],int nc,
					double UNITT,double k_B,double pi,
					double* EAlpha,double* EZ);

#endif
