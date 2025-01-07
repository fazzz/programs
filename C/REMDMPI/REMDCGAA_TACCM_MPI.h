

#ifndef INCLUDE_REMD_TAM
#define INCLUDE_REMD_TAM

#include <netcdf.h>
#include "netcdf_mineL.h"

#include "FFL.h"

int MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_Amber_AAFF(int myrank, int num_procs,int tag, MPI_Status* status,int numEX,
						    ////////////// AA ////////////////////////////
						    double *crdAA,double *velAA, double *mass, int numatom,
						    double zetaAA,double V_zetaAA, double QAA,
						    struct potential e, struct force f, double NfKTAA,
						    double *avePEAA, double *aveKEAA,double *aveTAA,
						    double *varPEAA, double *varKEAA,double *varTAA, 
						    struct my_netcdf_out_id_MCD nc_id_MCDAA,  FILE *outputfileAA, 
						    ////////////// CG ////////////////////////////
						    double *crdCG,double *velCG,
						    double zetaCG,double V_zetaCG, double QCG,
						    struct potential_GOLMAA_PROTEINS2008 e_GOLM, double NfKTCG,
						    double *avePECG, double *aveKECG,double *aveTCG,
						    double *varPECG, double *varKECG,double *varTCG, 
						    struct my_netcdf_out_id_MCD nc_id_MCDCG,  FILE *outputfileCG, 
						    ////////////// COMMON /////////////////////////
						    double dt,double dt2,double wdt2[3],double wdt4[3], int nc,
						    int numstep, int interval, double TAA, double TCG,
						    double UNITT, double k_B, double tau, double pi,
						    //////////////// TACCM ///////////////////////
						    double **Z,double **velZ,double massZ,
						    double *zetaZ,double *V_zetaZ,
						    double T0Z,double QZ,double NfKTZ,int numZ,
						    double Kapa,int **pairs,
						    double *avePEZ, double *aveKEZ,double *aveTZ,
						    double *varPEZ, double *varKEZ,double *varTZ, 
						    FILE **trjfileZ, FILE **trjfilTheta
						    //////////////// TACCM ///////////////////////
						    );

#endif
