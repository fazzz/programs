
#ifndef INCLUDE_REMD
#define INCLUDE_REMD

#include <netcdf.h>
#include "netcdf_mineL.h"

#include "FFL.h"

int MPI_TREMD_pep_NHC_MP1998_Amber_AAFF(int myrank, int num_procs,int tag, MPI_Status* status,
					int numEX,  int numRE,
					double **crd,double **vel, double *mass, int numatom,
					double *zeta,double *V_zeta, double *Q,
					struct potential *e, struct force *f,
					double *T,double *NfKT, int numstep,int interval,
					double dt,double dt2,double wdt2[3],double wdt4[3], int nc,
					double *avePE, double *aveKE,double *aveT,
					double *varPE, double *varKE,double *varT, 
					double UNITT, double k_B, double tau, double pi,
					struct my_netcdf_out_id_MCD *nc_id_MCD,  FILE **outputfile );

void readInputs(FILE *inputfile, double **crd, double **vel, int numatom, double *T0);

void readOutputs(FILE *outputfile, char *trjfilename[100], char *outputfilename[100]);

#endif
