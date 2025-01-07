
#ifndef INCLUDE_REMD_TAM2
#define INCLUDE_REMD_TAM2

#include "REMDCGAA_TACCM_MPI_2_Amber_PROTEINS2008.h"

struct TACCMDataforREMD_A_P2008_Z_asCA{  
  int numZ;
  double *Z,*velZ,massZ;
  double zetaZ,V_zetaZ,Q;

  double T,NfKT;
  double *avePEZ,*aveKEZ,*aveTZ;
  double *varPEZ,*varKEZ,*varTZ;

  double KZAA,KZCG;  int *index;

  struct my_netcdf_out_id_MCD trjfileZ,trjfilThetaAA,trjfilThetaCG;

  double **predict,**correct;

  double sZ,s_velZ,gziZ,gzi_velZ,predict_gziZ[5],correct_gziZ[5],predict_sZ[5],correct_sZ[5];
  double To,KEo;

};

double **MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_Amber_PROTEINS2008_Z_asCA(int myrank,int num_procs,
									int tag, MPI_Status* status,
									int numRE, int numEX, 
									double *KZAA, double *KZCG,
									struct AADataforREMD_Amber AAData,
									struct CGDataforREMD_PROTEINS2008 CGData,
									struct TACCMDataforREMD_A_P2008_Z_asCA ZData,
									struct AACGCommonDataforREMD_A_P2008 CData,
									double T0AA,double T0CG, double T0Z, 
									int numstep, int interval, 
									double dt,double dt2,
									double wdt2[3],double wdt4[3], int nc,
									double UNITT, double k_B, double tau, 
									double pi,
									FILE* logfile );

double  CE_TACCM_CGAA_Z_asCA(double *crdAA,double *crdCG,double *Z,
			     int numatom,int numZ,
			     double KZAA,double KZCG,int *pairs,double pi,
			     double *EA,double *ECG,double *EZ);

#endif
