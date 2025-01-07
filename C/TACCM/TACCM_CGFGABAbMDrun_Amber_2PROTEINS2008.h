
#ifndef INCLUDE_TA_MD
#define INCLUDE_TA_MD

#include <netcdf.h>
#include "netcdf_mineL.h"

#include "FFL.h"

#include "ABAb.h"

double runTACCM_2CG1FG_ABAbMD_NH_new_Amber_PROTEINS2008(// FG ////////////////////////////////////////////
							double *FGcrd,double *FGq,double *FGqvel,
							double *FGpredict,double *FGcorrect,
							double FGs, double FGs_vel, 
							double FGpredict_s[5],   double FGcorrect_s[5],
							double FGgzi, double FGgzi_vel,
							double FGpredict_gzi[5], double FGcorrect_gzi[5],
							int *FGnumclutparent,int *FGterminal,int *FGorigin,
							double *FGvel_Term,
							double **FGpredict_Term,double **FGpredict_Term2,
							double **FGcorrect_Term,double **FGcorrect_Term2,
							CLTb *FGclt, double FGQ,
							struct potential e, struct force f, 
							double TFG, double TFGo, double KEFGo,
							double *avePEFG, double *aveKEFG,double *aveTFG,
							double *varPEFG, double *varKEFG,double *varTFG,
							struct my_netcdf_out_id_MCD nc_id_MCDFG,  FILE *outputfileFG,
							// CG1 ////////////////////////////////////////////
							double *CG1crd,double *CG1q,double *CG1qvel,
							double *CG1predict,double *CG1correct,
							double CG1s, double CG1s_vel, 
							double CG1predict_s[5],   double CG1correct_s[5],
							double CG1gzi, double CG1gzi_vel,
							double CG1predict_gzi[5], double CG1correct_gzi[5],
							int *CG1numclutparent,int *CG1terminal,int *CG1origin,
							double *CG1vel_Term,
							double **CG1predict_Term,double **CG1predict_Term2,
							double **CG1correct_Term,double **CG1correct_Term2,
							CLTb *CG1clt, double CG1Q,
							struct potential_GOLMAA_PROTEINS2008 e_CG1, 
							double TCG1, double TCG1o, double KECG1o,
							double *avePECG1, double *aveKECG1,double *aveTCG1,
							double *varPECG1, double *varKECG1,double *varTCG1,
							struct my_netcdf_out_id_MCD nc_id_MCDCG1,  FILE *outputfileCG1,
							// CG2 ////////////////////////////////////////////
							double *CG2crd,double *CG2q,double *CG2qvel,
							double *CG2predict,double *CG2correct,
							double CG2s, double CG2s_vel, 
							double CG2predict_s[5],   double CG2correct_s[5],
							double CG2gzi, double CG2gzi_vel,
							double CG2predict_gzi[5], double CG2correct_gzi[5],
							int *CG2numclutparent,int *CG2terminal,int *CG2origin,
							double *CG2vel_Term,
							double **CG2predict_Term,double **CG2predict_Term2,
							double **CG2correct_Term,double **CG2correct_Term2,
							CLTb *CG2clt, double CG2Q,
							struct potential_GOLMAA_PROTEINS2008 e_CG2, 
							double TCG2, double TCG2o, double KECG2o,
							double *avePECG2, double *aveKECG2,double *aveTCG2,
							double *varPECG2, double *varKECG2,double *varTCG2,
							struct my_netcdf_out_id_MCD nc_id_MCDCG2,  FILE *outputfileCG2,
							// Z  /////////////////////////////////////////////
							double *Z,double *velZ,double massZ,
							double **predict_Z,double **correct_Z,
							double sZ,double s_velZ,double gziZ,double gzi_velZ,
							double predict_gziZ[5],double correct_gziZ[5],
							double predict_sZ[5],double correct_sZ[5],
							double QZ,double TZ, double TZo, double KEZo,
							int numZ,double KZFG, 
							double KZCG1, double KZCG1,int **pairs,
							double *avePEZ, double *aveKEZ,double *aveTZ,
							double *varPEZ, double *varKEZ,double *varTZ, 
							FILE *trjfileZ, FILE *trjfilThetaFG, 
							FILE *trjfilThetaCG1,FILE *trjfilThetaCG2,
							// CM  ///////////////////////////////////////////////
							double *mass, int numatom, int numclut,int DOF,
							int numstep, int interval,int *l,
							double dt,double tau,double tau2,
							double UNITT, double k_B,double pi,
							double *PEZFG, 
							double *PEZCG1, double *PEZCG2, double *PEZ);

#endif
