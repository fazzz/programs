
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "FFLc.h"
#include "PTLb.h"

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

#include "TACCM_CGAAMDrun_test_CGb.h"
#include "TACCM_MDrun.h"
#include "TACCM_MD.h"

#include "MDb_NHC_MP1996.h"
#include "MDrunb.h"
#include "MD.h"

#include "MB.h"
#include "LA.h"
#include "TOPO.h"
#include "mymath.h"

double ffLc_calcffandforce_AA(double *crd, int numatom,struct potential *ene,struct force *f,struct AmberParmL ap) {
  int i;
  int numnb,num14;
  double *n_d;

  numnb=(*ene).parm.numnb;
  num14=(*ene).parm.num14;

  ffLc_calcDIHE((*ene).p_d,n_d,crd,1,0,0,ap);
  ffLc_calcDIHE_force_Cartesian((*f).f_d,crd,ap); // 1111

  ffLc_calcANGLE((*ene).p_a,crd,ap);
  ffLc_calcANGLE_force_Cartesian((*f).f_a,crd,ap);

  ffLc_calcBOND((*ene).p_b,crd,ap);
  ffLc_calcBOND_force_Cartesian((*f).f_b,crd,ap);
  
  (*ene).p_t=0.0;
  (*ene).p_e_t=0.0;
  (*ene).p_LJ_t=0.0;
  (*ene).p_e_14_t=0.0;
  (*ene).p_LJ_14_t=0.0;
  (*ene).p_d_t=0.0;
  (*ene).p_a_t=0.0;
  (*ene).p_b_t=0.0;
  for (i=0;i<numatom*3;++i) (*f).f_t[i]=0.0;

  for (i=0;i<ap.NPHIH+ap.MPHIA;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }
  for (i=0;i<ap.NTHETH+ap.MTHETA;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];             // 0911
  }                                          // 0911
  for (i=0;i<ap.NBONH+ap.MBONA;++i) {        // 0911
    (*ene).p_t+=(*ene).p_b[i];               // 0911
    (*ene).p_b_t+=(*ene).p_b[i];
  }
  for (i=0;i<numatom*3;++i) {
    (*f).f_e_14[i]=1.0/1.2*(*f).f_e_14[i];
    (*f).f_LJ_14[i]=0.5*(*f).f_LJ_14[i];
  }

  for (i=0;i<numatom*3;++i)
    (*f).f_t[i]+=(*f).f_d[i]+(*f).f_a[i]+(*f).f_b[i];
  
  return (*ene).p_t;
}


