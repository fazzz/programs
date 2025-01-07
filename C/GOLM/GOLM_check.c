#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLM.h"
#include "GOLM_check.h"

#define ON 1
#define OFF 0

double *GOLMff_calcff_check(double *crd, int numatom,struct potential_GOLM *ene,
			   int flagb,int flaga, int flagd, int flagnc, int flagnn,
			   int numspatom,double dx) {
  int i,j;
  double *crddx,*crddy,*crddz;
  double *p_b_dx,*p_b_dy,*p_b_dz;
  double *p_a_dx,*p_a_dy,*p_a_dz;
  double *p_d_dx,*p_d_dy,*p_d_dz;
  double *p_natatt_dx,*p_natatt_dy,*p_natatt_dz;
  double *p_repul_dx,*p_repul_dy,*p_repul_dz;
  double f_b[3],f_a[3],f_d[3],f_natatt[3],f_repul[3];
  double *f;

  crddx=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddy=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddz=(double *)gcemalloc(sizeof(double)*numatom*3);

  p_b_dx=(double *)gcemalloc(sizeof(double)*(*ene).num_bond);
  p_b_dy=(double *)gcemalloc(sizeof(double)*(*ene).num_bond);
  p_b_dz=(double *)gcemalloc(sizeof(double)*(*ene).num_bond);
  p_a_dx=(double *)gcemalloc(sizeof(double)*(*ene).num_angl);
  p_a_dy=(double *)gcemalloc(sizeof(double)*(*ene).num_angl);
  p_a_dz=(double *)gcemalloc(sizeof(double)*(*ene).num_angl);
  p_d_dx=(double *)gcemalloc(sizeof(double)*(*ene).num_dihe);
  p_d_dy=(double *)gcemalloc(sizeof(double)*(*ene).num_dihe);
  p_d_dz=(double *)gcemalloc(sizeof(double)*(*ene).num_dihe);
  p_natatt_dx=(double *)gcemalloc(sizeof(double)*numatom);
  p_natatt_dy=(double *)gcemalloc(sizeof(double)*numatom);
  p_natatt_dz=(double *)gcemalloc(sizeof(double)*numatom);
  p_repul_dx=(double *)gcemalloc(sizeof(double)*numatom);
  p_repul_dy=(double *)gcemalloc(sizeof(double)*numatom);
  p_repul_dz=(double *)gcemalloc(sizeof(double)*numatom);

  f=(double *)gcemalloc(sizeof(double)*3);

  /**********************************/
  /* memcpy(crddx,crd,sizeof(crd)); */
  /* memcpy(crddy,crd,sizeof(crd)); */
  /* memcpy(crddz,crd,sizeof(crd)); */
  /**********************************/

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      crddx[i*3+j]=crd[i*3+j];
      crddy[i*3+j]=crd[i*3+j];
      crddz[i*3+j]=crd[i*3+j];
    }
  }

  crddx[numspatom*3]+=dx;
  crddy[numspatom*3+1]+=dx;
  crddz[numspatom*3+2]+=dx;

  if (flagb==ON) {
    GOLMpote_calcBOND((*ene).p_b,crd,(*ene).FC_bond,(*ene).BEQ,(*ene).num_bond,(*ene).index_bond);
    GOLMpote_calcBOND(p_b_dx,crddx,(*ene).FC_bond,(*ene).BEQ,(*ene).num_bond,(*ene).index_bond);
    GOLMpote_calcBOND(p_b_dy,crddy,(*ene).FC_bond,(*ene).BEQ,(*ene).num_bond,(*ene).index_bond);
    GOLMpote_calcBOND(p_b_dz,crddz,(*ene).FC_bond,(*ene).BEQ,(*ene).num_bond,(*ene).index_bond);

    for (i=0;i<3;++i) {
      f_b[i]=0.0;
    }
    for (i=0;i<(*ene).num_bond;++i) {
      f_b[0]+=(p_b_dx[i]-(*ene).p_b[i])/dx*UNIT;
      f_b[1]+=(p_b_dy[i]-(*ene).p_b[i])/dx*UNIT;
      f_b[2]+=(p_b_dz[i]-(*ene).p_b[i])/dx*UNIT;
    }
  }

  if (flaga==ON) {
    GOLMpote_calcANGLE((*ene).p_a,crd,(*ene).FC_angle,(*ene).AEQ,(*ene).num_angl,(*ene).index_angl);
    GOLMpote_calcANGLE(p_a_dx,crddx,(*ene).FC_angle,(*ene).AEQ,(*ene).num_angl,(*ene).index_angl);
    GOLMpote_calcANGLE(p_a_dy,crddy,(*ene).FC_angle,(*ene).AEQ,(*ene).num_angl,(*ene).index_angl);
    GOLMpote_calcANGLE(p_a_dz,crddz,(*ene).FC_angle,(*ene).AEQ,(*ene).num_angl,(*ene).index_angl);

    for (i=0;i<3;++i) {
      f_a[i]=0.0;
    }
    for (i=0;i<(*ene).num_angl;++i) {
      f_a[0]+=(p_a_dx[i]-(*ene).p_a[i])/dx*UNIT;
      f_a[1]+=(p_a_dy[i]-(*ene).p_a[i])/dx*UNIT;
      f_a[2]+=(p_a_dz[i]-(*ene).p_a[i])/dx*UNIT;
    }
  }

  if (flagd==ON) {
    GOLMpote_calcDIHE((*ene).p_d,crd,(*ene).DEQ,(*ene).FC_dihed1,(*ene).FC_dihed2,(*ene).num_dihe,(*ene).index_dihe);
    GOLMpote_calcDIHE(p_d_dx,crddx,(*ene).DEQ,(*ene).FC_dihed1,(*ene).FC_dihed2,(*ene).num_dihe,(*ene).index_dihe);
    GOLMpote_calcDIHE(p_d_dy,crddy,(*ene).DEQ,(*ene).FC_dihed1,(*ene).FC_dihed2,(*ene).num_dihe,(*ene).index_dihe);
    GOLMpote_calcDIHE(p_d_dz,crddz,(*ene).DEQ,(*ene).FC_dihed1,(*ene).FC_dihed2,(*ene).num_dihe,(*ene).index_dihe);

    for (i=0;i<3;++i) {
      f_d[i]=0.0;
    }
    for (i=0;i<(*ene).num_dihe;++i) {
      f_d[0]+=(p_d_dx[i]-(*ene).p_d[i])/dx*UNIT;
      f_d[1]+=(p_d_dy[i]-(*ene).p_d[i])/dx*UNIT;
      f_d[2]+=(p_d_dz[i]-(*ene).p_d[i])/dx*UNIT;
    }
  }

  if (flagnc==ON) {
    GOLMpote_calcNatAtt((*ene).p_natatt,crd,(*ene).ALJ_natatt,(*ene).BLJ_natatt,(*ene).ep_natatt,(*ene).num_natatt,numatom,(*ene).index_natatt);
    GOLMpote_calcNatAtt(p_natatt_dx,crddx,(*ene).ALJ_natatt,(*ene).BLJ_natatt,(*ene).ep_natatt,(*ene).num_natatt,numatom,(*ene).index_natatt);
    GOLMpote_calcNatAtt(p_natatt_dy,crddy,(*ene).ALJ_natatt,(*ene).BLJ_natatt,(*ene).ep_natatt,(*ene).num_natatt,numatom,(*ene).index_natatt);
    GOLMpote_calcNatAtt(p_natatt_dz,crddz,(*ene).ALJ_natatt,(*ene).BLJ_natatt,(*ene).ep_natatt,(*ene).num_natatt,numatom,(*ene).index_natatt);

    f_natatt[0]=(p_natatt_dx[numspatom]-(*ene).p_natatt[numspatom])/dx*UNIT;
    f_natatt[1]=(p_natatt_dy[numspatom]-(*ene).p_natatt[numspatom])/dx*UNIT;
    f_natatt[2]=(p_natatt_dz[numspatom]-(*ene).p_natatt[numspatom])/dx*UNIT;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////

  if (flagnn==ON) {
    GOLMpote_calcRepul((*ene).p_repul,crd,(*ene).ALJ_repul,numatom);
    GOLMpote_calcRepul(p_repul_dx,crddx,(*ene).ALJ_repul,numatom);
    GOLMpote_calcRepul(p_repul_dy,crddy,(*ene).ALJ_repul,numatom);
    GOLMpote_calcRepul(p_repul_dz,crddz,(*ene).ALJ_repul,numatom);

    f_repul[0]=(p_repul_dx[numspatom]-(*ene).p_repul[numspatom])/dx*UNIT;
    f_repul[1]=(p_repul_dy[numspatom]-(*ene).p_repul[numspatom])/dx*UNIT;
    f_repul[2]=(p_repul_dz[numspatom]-(*ene).p_repul[numspatom])/dx*UNIT;
  }

  f[0]=0.0;
  f[1]=0.0;
  f[2]=0.0;
  if(flagb==ON) {
    f[0]+=f_b[0];
    f[1]+=f_b[1];
    f[2]+=f_b[2];
  }
  if(flaga==ON) {
    f[0]+=f_a[0];
    f[1]+=f_a[1];
    f[2]+=f_a[2];
  }
  if (flagd==ON) {
    f[0]+=f_d[0];
    f[1]+=f_d[1];
    f[2]+=f_d[2];
  }
  if (flagnc==ON) {
    f[0]+=f_natatt[0];
    f[1]+=f_natatt[1];
    f[2]+=f_natatt[2];
  }
  if (flagnn==ON) {
    f[0]+=f_repul[0];
    f[1]+=f_repul[1];
    f[2]+=f_repul[2];
  }

  return f;
}
