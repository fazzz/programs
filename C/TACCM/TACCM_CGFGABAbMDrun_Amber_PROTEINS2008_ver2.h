
#ifndef INCLUDE_TA_MD
#define INCLUDE_TA_MD

#include <netcdf.h>
#include "netcdf_mineL.h"

#include "FFL.h"

#include "ABAb.h"

double runTACCM_CGFG_ABAbMD_NH_new_Amber_PROTEINS2008_ver2(// FG ////////////////////////////////////////////
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
							   struct my_netcdf_out_id_MCD nc_id_MCDFG,
							   FILE *outputfileFG,
							   // CG ////////////////////////////////////////////
							   double *CGcrd,double *CGq,double *CGqvel,
							   double *CGpredict,double *CGcorrect,
							   double CGs, double CGs_vel, 
							   double CGpredict_s[5],   double CGcorrect_s[5],
							   double CGgzi, double CGgzi_vel,
							   double CGpredict_gzi[5], double CGcorrect_gzi[5],
							   int *CGnumclutparent,int *CGterminal,int *CGorigin,
							   double *CGvel_Term,
							   double **CGpredict_Term,double **CGpredict_Term2,
							   double **CGcorrect_Term,double **CGcorrect_Term2,
							   CLTb *CGclt, double CGQ,
							   struct potential_GOLMAA_PROTEINS2008 e_CG, 
							   double TCG, double TCGo, double KECGo,
							   double *avePECG, double *aveKECG,double *aveTCG,
							   double *varPECG, double *varKECG,double *varTCG,
							   struct my_netcdf_out_id_MCD nc_id_MCDCG,
							   FILE *outputfileCG,
							   // Z  /////////////////////////////////////////////
							   double *Zcrd,double *Zq,double *Zqvel,
							   double *Zpredict,double *Zcorrect,
							   double Zs, double Zs_vel, 
							   double Zpredict_s[5],   double Zcorrect_s[5],
							   double Zgzi, double Zgzi_vel,
							   double Zpredict_gzi[5], double Zcorrect_gzi[5],
							   int *Znumclutparent,int *Zterminal,int *Zorigin,
							   double *Zvel_Term,
							   double **Zpredict_Term,double **Zpredict_Term2,
							   double **Zcorrect_Term,double **Zcorrect_Term2,
							   CLTb *Zclt, double ZQ,
							   double TZ, double TZo, double KEZo,
							   double *Z, int numZ, double massZ,
							   double KZFG, double KZCG,int **pairs,
							   double *avePEZ, double *aveKEZ,double *aveTZ,
							   double *varPEZ, double *varKEZ,double *varTZ, 
							   struct my_netcdf_out_id_MCD nc_id_MCDZ,
							   FILE *trjfileZ, FILE *trjfilThetaFG, FILE *trjfilThetaCG,
							   // CM  ///////////////////////////////////////////////
							   double *mass, int numatom, int numclut,int DOF,
							   int numstep, int interval,int *l,
							   double dt,double tau,double tau2,
							   double UNITT, double k_B,double pi,
							   double *PEZFG, double *PEZCG, double *PEZ);

#endif
