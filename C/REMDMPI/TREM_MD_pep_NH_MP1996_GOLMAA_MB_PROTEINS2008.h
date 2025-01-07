#ifndef INCLUDE_TREMD_GOLMAA_MB
#define INCLUDE_TREMD_GOLMAA_MB

#include <netcdf.h>
#include "netcdf_mineL.h"

#include "FFL.h"

int MPI_TREMD_pep_NH_MP1996_GOLMAA_MB_PROTEINS2008(int myrank, int num_procs,int tag, MPI_Status* status,
						   int numEX,  int numRE,
						   double **crd,double **vel, double *mass, 
						   int numatom, int numheavyatom,
						   double *zeta,double *V_zeta, double *Q,
						   struct potential_GOLMAA_MB_PROTEINS2008 *e_GOLM,
						   struct potential *e,
						   double de, double d2,
						   double *T,double *NfKT, int numstep,int interval,
						   double dt,double dt2,double wdt2[3],double wdt4[3], int nc,
						   double *avePE, double *aveKE,double *aveT,
						   double *varPE, double *varKE,double *varT, 
						   double UNITT, double k_B, double tau, double pi,
						   struct my_netcdf_out_id_MCD *nc_id_MCD,
						   FILE **outputfile, FILE *logfile );

void readInputs_TREM_GOLMAA(FILE *inputfile, double **crd, double **vel, int numatom, double *T0, int numRE);

void readOutputs_readInputs_TREM_GOLMAA(FILE *outputfile, char *trjfilename[100], char *outputfilename[100]);

double runMD_pep_NH_MP1996_GOLMAA_MB_PROTEINS2008(double *crd,double *vel, double *mass, 
						  int numatom, int numheavyatom,
						  double *zeta,double *V_zeta, double Q,
						  struct potential_GOLMAA_MB_PROTEINS2008 e_GOLM, 
						  struct potential e,
						  double de, double d2,
						  double T, double NfKT, int numstep, int interval,int *l,
						  double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
						  double *avePE, double *aveKE,double *aveT,
						  double *varPE, double *varKE,double *varT,
						  double UNITT, double k_B,
						  struct my_netcdf_out_id_MCD nc_id_MCD,  FILE *outputfile);

#endif


