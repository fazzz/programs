

#ifndef INCLUDE_REMD_TAM
#define INCLUDE_REMD_TAM

#include <netcdf.h>
#include "mpi.h"

#include "REMDMPI.h"
#include "REMD_TACCM_MPI.h"

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

#include "FFL.h"
#include "PTL.h"
#include "EF.h"

#include "RAND.h"
#include "BOXMULL.h"
#include "MC.h"

#include "MDrun.h"
#include "TACCM.h"
#include "TACCM_CGAAMDrun_Amber_PROTEINS2008.h"
//#include "TACCM_CGAAMDrun_test_CG.h"

#include "netcdf_mineL.h"

#include "ABAb.h"

#define AAINPF 0
#define CGINPF 1
#define AAKZ 2
#define CGKZ 3

struct AADataforREMD_Amber{  
  double *crd,*vel;
  double zeta,V_zeta,Q;
  struct potential e;
  struct force f;
  double T,NfKT;
  double *avePE,*aveKE,*aveT;
  double *varPE,*varKE,*varT;

  CLTb *clt;
  double *q,*qvel;
  double *predict,*correct;
  double s,s_vel,s_acc,gzi,gzi_vel,predict_gzi[5],correct_gzi[5],predict_s[5],correct_s[5];
  double *vel_Term,**predict_Term,**predict_Term2,**correct_Term,**correct_Term2;
  double To,KEo;

  struct my_netcdf_out_id_MCD nc_id_MCD;  FILE *outputfile;
};

struct CGDataforREMD_PROTEINS2008{  
  double *crd,*vel;
  double zeta,V_zeta,Q;
  struct potential_GOLMAA_PROTEINS2008 e;
  double T,NfKT;
  double *avePE,*aveKE,*aveT;
  double *varPE,*varKE,*varT;

  CLTb *clt;
  double *q,*qvel;
  double *predict,*correct;
  double s,s_vel,s_acc,gzi,gzi_vel,predict_gzi[5],correct_gzi[5],predict_s[5],correct_s[5];
  double *vel_Term,**predict_Term,**predict_Term2,**correct_Term,**correct_Term2;
  double To,KEo;

  struct my_netcdf_out_id_MCD nc_id_MCD;  FILE *outputfile;
};

struct AACGCommonDataforREMD_A_P2008{  
  int numatom,numheavyatom,numres;

  double *mass;  int numstep,interval;

  int numclut;

  int *numclutparent,*terminal,*origin;
  int DOF;

};

struct TACCMDataforREMD_A_P2008{  
  int numZ;
  double *Z,*velZ,massZ;
  double zetaZ,V_zetaZ;

  double T,QZ,NfKTZ;
  double KZAA,KZCG;  int **pairs;

  double *avePEZ,*aveKEZ,*aveTZ;
  double *varPEZ,*varKEZ,*varTZ;

  FILE *trjfileZ,*trjfilThetaAA,*trjfilThetaCG;

  double **predict,**correct;

  double sZ,s_velZ,gziZ,gzi_velZ,predict_gziZ[5],correct_gziZ[5],predict_sZ[5],correct_sZ[5];
  double To,KEo;

};

struct TACCMDataforREMD_A_2P2008{  
  int numZ;
  double *Z,*velZ,massZ;
  double zetaZ,V_zetaZ;

  double T,QZ,NfKTZ;
  double KZAA,KZCG[2];  int **pairs;

  double *avePEZ,*aveKEZ,*aveTZ;
  double *varPEZ,*varKEZ,*varTZ;

  FILE *trjfileZ,*trjfilThetaAA,*trjfilThetaCG1,*trjfilThetaCG2;

  double **predict,**correct;

  double sZ,s_velZ,gziZ,gzi_velZ,predict_gziZ[5],correct_gziZ[5],predict_sZ[5],correct_sZ[5];
  double To,KEo;

};

struct TACCMDataforREMD_A_P2008_ver2{  
  int numZ;
  double *Z,massZ;

  double KZAA,KZCG;
  int **pairs;

  double *crd,*vel;
  double zeta,V_zeta,Q;
  struct potential e;
  struct force f;
  double T,NfKT;
  double *avePE,*aveKE,*aveT;
  double *varPE,*varKE,*varT;

  CLTb *clt;
  double *q,*qvel;
  double *predict,*correct;
  double s,s_vel,s_acc,gzi,gzi_vel,predict_gzi[5],correct_gzi[5],predict_s[5],correct_s[5];
  double *vel_Term,**predict_Term,**predict_Term2,**correct_Term,**correct_Term2;
  double To,KEo;

  struct my_netcdf_out_id_MCD nc_id_MCD;

  FILE *trjfileZ,*trjfilThetaAA,*trjfilThetaCG;

};

double **MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_Amber_PROTEINS2008(int myrank,int num_procs,
								 int tag, MPI_Status* status,
								 int numRE, int numEX, double *KZAA, double *KZCG,
								 struct AADataforREMD_Amber AAData,
								 struct CGDataforREMD_PROTEINS2008 CGData,
								 struct TACCMDataforREMD_A_P2008 ZData,
								 struct AACGCommonDataforREMD_A_P2008 CData,
								 double T0AA,double T0CG, double T0Z, 
								 int numstep, int interval, 
								 double dt,double dt2,
								 double wdt2[3],double wdt4[3], int nc,
								 double UNITT, double k_B, double tau, double pi,
								 FILE* logfile );

double **MPI_CGFGTREM_TACCM_ABAbMD_NH_new_Amber_PROTEINS2008(int myrank,int num_procs,
							     int tag, MPI_Status* status,
							     int numRE, int numEX, double *KZAA, double *KZCG,
							     struct AADataforREMD_Amber AAData,
							     struct CGDataforREMD_PROTEINS2008 CGData,
							     struct TACCMDataforREMD_A_P2008 ZData,
							     struct AACGCommonDataforREMD_A_P2008 CData,
							     double T0AA,double T0CG, double T0Z, 
							     int numstep, int interval, 
							     double dt,double tau, double tau2,
							     double UNITT, double k_B, double pi, FILE* logfile );

double **MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_2PROTEINS2008(int myrank,int num_procs,
							    int tag, MPI_Status* status,
							    int numRE, int numEX, double *KZCG1, double *KZCG2,
							    struct CGDataforREMD_PROTEINS2008 CG1Data,
							    struct CGDataforREMD_PROTEINS2008 CG2Data,
							    struct TACCMDataforREMD_A_P2008 ZData,
							    struct AACGCommonDataforREMD_A_P2008 CData,
							    double T0CG1,double T0CG2, double T0Z, 
							    int numstep, int interval, 
							    double dt,double dt2,
							    double wdt2[3],double wdt4[3], int nc,
							    double UNITT, double k_B, double tau, double pi,
							    FILE* logfile );

void  CGAAREMDreadInputs_test(FILE *inputfile,int numatom,int numRE,int myrank,
			      double *crdAA,double *velAA, double *crdCG,double *velCG,
			      double *KZAA, double *KZCG);

#endif
