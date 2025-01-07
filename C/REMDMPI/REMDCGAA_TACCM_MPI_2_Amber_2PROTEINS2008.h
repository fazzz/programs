

#ifndef INCLUDE_REMD_TAM2
#define INCLUDE_REMD_TAM2

#include "REMDCGAA_TACCM_MPI_2_Amber_PROTEINS2008.h"

#define AA1INPF 0
#define CG1INPF 1
#define CG2INPF 2
#define AA1KZ 3
#define CG1KZ 4
#define CG2KZ 5

double **MPI_CGFGTREM_TACCM_ABAbMD_NH_new_Amber_2PROTEINS2008(int myrank,int num_procs,
							      int tag, MPI_Status* status,
							      int numRE, int numEX, double *KZAA, double **KZCG,
							      struct AADataforREMD_Amber AAData,
							      struct CGDataforREMD_PROTEINS2008 CGData[2],
							      struct TACCMDataforREMD_A_2P2008 ZData,
							      struct AACGCommonDataforREMD_A_P2008 CData,
							      double T0AA,double T0CG[2], double T0Z, 
							      int numstep, int interval, 
							      double dt,double tau, double tau2,
							      double UNITT, double k_B, double pi, FILE* logfile );

void  CGAAREMDreadInputs_Amber_PROTEINS2008_Amber_hybrid_1FG2CG(FILE *inputfile,int numatom,int numRE,int myrank,
								double *crdAA,double *velAA, 
								double *crdCG1,double *velCG1,
								double *crdCG2,double *velCG2,
								double *KZAA, double *KZCG1, double *KZCG2);

double CE_TACCM_2CG1AA(double *crdAA,double *crdCG1,double *crdCG2,double *Z, int numatom,int numZ,
		       double KZAA,double KZCG1,double KZCG2,int **pairs,double pi,
		       double *EAA,double *ECG1,double *ECG2,double *EZ);

double **MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_Amber_2PROTEINS2008(int myrank,int num_procs,
								  int tag, MPI_Status* status,
								  int numRE, int numEX, double *KZAA, double **KZCG,
								  struct AADataforREMD_Amber AAData,
								  struct CGDataforREMD_PROTEINS2008 CGData[2],
								  struct TACCMDataforREMD_A_2P2008 ZData,
								  struct AACGCommonDataforREMD_A_P2008 CData,
								  double T0AA,double T0CG[2], double T0Z, 
								  int numstep, int interval, 
								  double dt,double dt2,
								  double wdt2[3],double wdt4[3], int nc,
								  double UNITT, double k_B, double tau, double pi,
								  FILE* logfile );

#endif
