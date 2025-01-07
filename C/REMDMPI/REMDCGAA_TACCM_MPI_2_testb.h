

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

#define AAINPF 0
#define CGINPF 1
#define AAKZ 2
#define CGKZ 3

struct AADataforREMD_testb{  
  double *crd,*vel;
  double zeta,V_zeta,Q;
  struct potential e;
  struct force f;
  double T,NfKT;
  double *avePE,*aveKE,*aveT;
  double *varPE,*varKE,*varT;

  struct my_netcdf_out_id_MCD nc_id_MCD;  FILE *outputfile;
};

struct CGDataforREMD_testb{  
  double *crd,*vel;
  double zeta,V_zeta,Q;
  struct potential_GOLMAA_PROTEINS2008 e_GOLM;
  double T,NfKT;
  double *avePE,*aveKE,*aveT;
  double *varPE,*varKE,*varT;

  struct my_netcdf_out_id_MCD nc_id_MCD;  FILE *outputfile;
};

struct AACGCommonDataforREMD_testb{  
  int numatom,numheavyatom,numres;

  double *mass;  int numstep,interval;
};

struct TACCMDataforREMD_testb{  
  int numZ;
  double *Z,*velZ,massZ;
  double zetaZ,V_zetaZ;

  double T,QZ,NfKTZ;
  double KZAA,KZCG;  int **pairs;

  double *avePEZ,*aveKEZ,*aveTZ;
  double *varPEZ,*varKEZ,*varTZ;

  FILE *trjfileZ,*trjfilThetaAA,*trjfilThetaCG;
};

double **MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_testb(int myrank,int num_procs,int tag, MPI_Status* status,
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

void  CGAAREMDreadInputs_testb(FILE *inputfile,int numatom,int numRE,int myrank,
			      double *crdAA,double *velAA, double *crdCG,double *velCG,
			      double *KZAA, double *KZCG);

#endif
