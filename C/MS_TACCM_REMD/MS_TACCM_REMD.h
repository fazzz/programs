
#ifndef INCLUDE_REMD_TAM
#define INCLUDE_REMD_TAM

#include <netcdf.h>
#include "netcdf_mineL.h"

#include "FFLc.h"
#include "PTLb.h"

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

struct AADataMSREMD_Amber{  
  double *crd,*vel;
  double zeta,V_zeta,Q;
  struct potential e;
  struct force f;
  double T,NfKT;
  double *avePE,*aveKE,*aveT;
  double *varPE,*varKE,*varT;

  struct my_netcdf_out_id_MCD nc_id_MCD;  FILE *outputfile;
};

struct CGDataMSREMD_PROTEINS2008{  
  double *crd,*vel;
  double zeta,V_zeta,Q;
  struct potential_GOLMAA_PROTEINS2008 e_GOLM;
  double T,NfKT;
  double *avePE,*aveKE,*aveT;
  double *varPE,*varKE,*varT;

  struct my_netcdf_out_id_MCD nc_id_MCD;  FILE *outputfile;
};

struct ZDataMSREMD_A_P2008_pca{  
  int numZ;
  int numatom;
  double *Z,*velZ,massZ;
  double zetaZ,V_zetaZ;

  double T,QZ,NfKTZ;
  double KZAA,KZCG;

  double *A;

  double *mean;
  double *v,*E;

  double *avePEZ,*aveKEZ,*aveTZ;
  double *varPEZ,*varKEZ,*varTZ;

  FILE *trjfileZ;
};

double run_multi_MS_TACCM_Amber_PROTEINS2008_pca(int numAA, int numCG,
						 struct AADataMSREMD_Amber *FGdata,
						 struct CGDataMSREMD_PROTEINS2008 *CGdata,
						 struct ZDataMSREMD_A_P2008_pca Zdata,
						 struct AmberParmL *ap,struct AmberParmL *ap_CG,
						 double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
						 double UNITT, double k_B,double pi,
						 double *PEZAA, double *PEZCG,double *PEZ);


#endif
