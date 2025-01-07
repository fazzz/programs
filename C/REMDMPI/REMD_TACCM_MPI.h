
#ifndef INCLUDE_REMD_TAM
#define INCLUDE_REMD_TAM

#include <netcdf.h>
#include "netcdf_mineL.h"

#include "FFL.h"

int MPI_TREM_TACCM_MD_pep_NHC_MP1998_Amber_AAFF(int myrank, int num_procs,int tag, MPI_Status* status,
						int numEX,  int numRE,
						double **crd,double **vel, double *mass, int numatom,
						double *zeta,double *V_zeta, double Q,
						struct potential *e, struct force *f,
						double T,double NfKT, int numstep,int interval,
						double dt,double dt2,double wdt2[3],double wdt4[3], int nc,
						double *avePE, double *aveKE,double *aveT,
						double *varPE, double *varKE,double *varT, 
						double UNITT, double k_B, double tau, double pi,
						struct my_netcdf_out_id_MCD *nc_id_MCD,  FILE **outputfile, 
						//////////////// TACCM ///////////////////////
						double **Z,double **velZ,double massZ,
						double *zetaZ,double *V_zetaZ,
						double *TZ,double *QZ,double *NfKTZ,int numZ,
						double Kapa,int **pairs,
						double *avePEZ, double *aveKEZ,double *aveTZ,
						double *varPEZ, double *varKEZ,double *varTZ, 
						FILE **trjfileZ, FILE **trjfilTheta
						//////////////// TACCM ///////////////////////
						);

#endif
