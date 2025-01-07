
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLMAA_hybrid_set.h"
#include "GOLMAA_hybrid.h"

#include "PTL.h"

double GOLMAA_hyb_ff_calcff(double *crd, int numatom,struct potential_GOLMAA_hybrid *ene) {
  int i,j;

  GOLMAA_hyb_ff_calcff_nonlocal(crd,numatom,
				(*ene).NC_index,(*ene).numNC,(*ene).NotNC_index,(*ene).numNotNC,
				(*ene).ALJ_natatt,(*ene).BLJ_natatt,(*ene).ep_natatt,(*ene).ALJ_repul,
				(*ene).p_natatt,(*ene).p_repul,(*ene).f_natatt,(*ene).f_repul);

  (*ene).p_t=0.0;
  (*ene).p_natatt_t=0.0;
  (*ene).p_repul_t=0.0;
  for (i=0;i<numatom;++i)
    (*ene).p_natatt_t+=(*ene).p_natatt[i];
  for (i=0;i<numatom;++i)
    (*ene).p_repul_t+=(*ene).p_repul[i];
  (*ene).p_t=(*ene).p_natatt_t+(*ene).p_repul_t;

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) (*ene).f_t[i][j]=(*ene).f_natatt[i][j]+(*ene).f_repul[i][j];

}

double GOLMAA_hyb_ff_calcff_nonlocal(double *crd, int numatom,int *NC_index,int numNC,int *NotNC_index,int numNotNC,
				     double *ALJ_natatt,double *BLJ_natatt,double e_natatt,double ALJ_Repul,
				     double *p_natatt,double *p_repul,double **f_natatt,double **f_repul) {
  int i,j,k;
  int atomi,atomj;

  double vec[3];
  double len,len2,len6,len12;
  double p12,p6;

  for (i=0;i<numatom;++i) {
    p_natatt[i] = 0.0;
    p_repul[i] = 0.0;
  }

  for (i=0;i<numatom;++i) {
    for (k=0;k<3;++k) {
      f_natatt[i][k] = 0.0;
      f_repul[i][k] = 0.0;
    }
  }

  for (i=0;i<numNC;++i) {
    atomi=NC_index[i*2+0];
    atomj=NC_index[i*2+1];

    len2 = 0.0;
    for(k=0;k<3;++k){
      vec[k] = crd[atomi*3+k]-crd[atomj*3+k];
      len2 += vec[k]*vec[k];
    }
    len = sqrt(len2);
    len6=len2;
    len12=len2;
    for (k=0;k<2;++k)  len6 = len6*len2;
    for (k=0;k<5;++k)  len12 = len12*len2;
    p12 = ALJ_natatt[i]/len12;
    p6 = BLJ_natatt[i]/len6;
    p_natatt[atomi] += e_natatt*(p12-2.0*p6);
    //    p_natatt[j] = e_natatt*(p12-2.0*p6);
    for (k=0;k<3;++k) {
      f_natatt[atomi][k] +=  e_natatt*(12.0*p12-6.0*2.0*p6)/(len2)*vec[k]*UNIT;
      f_natatt[atomj][k] += -e_natatt*(12.0*p12-6.0*2.0*p6)/(len2)*vec[k]*UNIT;
    }
  }

  for (i=0;i<numNotNC;++i) {
    atomi=NotNC_index[i*2+0];
    atomj=NotNC_index[i*2+1];

    len2 = 0.0;
    for(k=0;k<3;++k){
      vec[k] = crd[atomi*3+k]-crd[atomj*3+k];
      len2 += vec[k]*vec[k];
    }
    len = sqrt(len2);
    len12=len2;
    for (k=0;k<5;++k) len12 = len12*len2;
    p12 = ALJ_Repul/len12;
    p_repul[atomi] += p12;
    //    p_repul[j] = p12;
    for (k=0;k<3;++k) {
      f_repul[atomi][k] +=  12.0*p12/(len2)*vec[k]*UNIT;
      f_repul[atomj][k] += -12.0*p12/(len2)*vec[k]*UNIT;
    }
  }
}

