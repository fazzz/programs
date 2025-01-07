
#ifndef INCLUDE_REMD_TAM_UENE
#define INCLUDE_REMD_TAM_UENE

#include "netcdf_mineL.h"

void  CGAAREMDreadInputs_calc_uene(FILE *inputfile,int numatom,int numRE,
				   double *KZAAs,double *KZCGs,
				   char **trjfilenameAA, char **trjfilenameCG,
				   struct my_netcdf_out_id_MCD *nc_id_MD_AA,
				   struct my_netcdf_out_id_MCD *nc_id_MD_CG,
				   FILE **trjfileZ);

void  CGAAREMDreadInputs_pickup_trj(FILE *inputfile, int numRE,
				    char **trjfilename, FILE **trjfile );

void  CGAAREMDreadInputs_calc_uene_1FG2CG(FILE *inputfile,int numatom,int numRE,
					  double *KZAAs,double *KZCG1s,double *KZCG2s,
					  char **trjfilenameAA, char **trjfilenameCG1, char **trjfilenameCG2,
					  struct my_netcdf_out_id_MCD *nc_id_MD_AA,
					  struct my_netcdf_out_id_MCD *nc_id_MD_CG1,
					  struct my_netcdf_out_id_MCD *nc_id_MD_CG2,
					  FILE **trjfileZ);

void  CGAAREMDreadInputs_calc_uene_InV2InW(FILE *inputfile,int numatom,int numRE,
					   double *KZAAs,double *KZCGs,
					   char **trjfilenameAA, char **trjfilenameCG,
					   struct my_netcdf_out_id_MCD *nc_id_MD_AA,
					   struct my_netcdf_out_id_MCD *nc_id_MD_CG,
					   FILE **trjfileZ,FILE **inputwfile);

void  CGAAREMDreadInputs_calc_uene_1FG2CG_InV2InW(FILE *inputfile,int numatom,int numRE,
						  double *KZAAs,double *KZCG1s,double *KZCG2s,
						  char **trjfilenameAA, char **trjfilenameCG1, char **trjfilenameCG2,
						  struct my_netcdf_out_id_MCD *nc_id_MD_AA,
						  struct my_netcdf_out_id_MCD *nc_id_MD_CG1,
						  struct my_netcdf_out_id_MCD *nc_id_MD_CG2,
						  FILE **trjfileZ,FILE **inputwfile);

#endif
