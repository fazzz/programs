#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLMAA.h"
#include "GOLMAA_check.h"
#include "GOLMAA_dbasin.h"

#define ON 1
#define OFF 0

double *GOLMAAff_calcff_check(double *crd, int numatom,struct potential_GOLMAA *ene,
			      int numspatom,double dx,double f_natatt[3],double f_repul[3],int **nb_matrix) {
  int i,j;
  double *crddx,*crddy,*crddz;
  double *p_natatt_dx,*p_natatt_dy,*p_natatt_dz;
  double *p_repul_dx,*p_repul_dy,*p_repul_dz;
  double *f;

  crddx=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddy=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddz=(double *)gcemalloc(sizeof(double)*numatom*3);

  p_natatt_dx=(double *)gcemalloc(sizeof(double)*numatom);
  p_natatt_dy=(double *)gcemalloc(sizeof(double)*numatom);
  p_natatt_dz=(double *)gcemalloc(sizeof(double)*numatom);
  p_repul_dx=(double *)gcemalloc(sizeof(double)*numatom);
  p_repul_dy=(double *)gcemalloc(sizeof(double)*numatom);
  p_repul_dz=(double *)gcemalloc(sizeof(double)*numatom);

  f=(double *)gcemalloc(sizeof(double)*3);

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

  GOLMAApoteforc_calcNatAttandRepul(p_natatt_dx,p_repul_dx,(*ene).f_natatt,(*ene).f_repul,crddx,(*ene).ALJ_natatt,(*ene).BLJ_natatt,(*ene).ep_natatt,(*ene).ALJ_repul,(*ene).num_natatt,numatom,(*ene).ncmap,nb_matrix);
  GOLMAApoteforc_calcNatAttandRepul(p_natatt_dy,p_repul_dy,(*ene).f_natatt,(*ene).f_repul,crddy,(*ene).ALJ_natatt,(*ene).BLJ_natatt,(*ene).ep_natatt,(*ene).ALJ_repul,(*ene).num_natatt,numatom,(*ene).ncmap,nb_matrix);
  GOLMAApoteforc_calcNatAttandRepul(p_natatt_dz,p_repul_dz,(*ene).f_natatt,(*ene).f_repul,crddz,(*ene).ALJ_natatt,(*ene).BLJ_natatt,(*ene).ep_natatt,(*ene).ALJ_repul,(*ene).num_natatt,numatom,(*ene).ncmap,nb_matrix);

  f_natatt[0]=-(p_natatt_dx[numspatom]-(*ene).p_natatt[numspatom])/dx*UNIT;
  f_natatt[1]=-(p_natatt_dy[numspatom]-(*ene).p_natatt[numspatom])/dx*UNIT;
  f_natatt[2]=-(p_natatt_dz[numspatom]-(*ene).p_natatt[numspatom])/dx*UNIT;

  f_repul[0]=-(p_repul_dx[numspatom]-(*ene).p_repul[numspatom])/dx*UNIT;
  f_repul[1]=-(p_repul_dy[numspatom]-(*ene).p_repul[numspatom])/dx*UNIT;
  f_repul[2]=-(p_repul_dz[numspatom]-(*ene).p_repul[numspatom])/dx*UNIT;

  f[0]=0.0;
  f[1]=0.0;
  f[2]=0.0;
  f[0]+=f_natatt[0];
  f[1]+=f_natatt[1];
  f[2]+=f_natatt[2];
  f[0]+=f_repul[0];
  f[1]+=f_repul[1];
  f[2]+=f_repul[2];

  return f;
}

double *GOLMAAff_dbasin_calcff_check(double *crd, int numatom,int numspatom,double dx,
				     double *refcrd1,double *refcrd2,
				     double R_C_D,double delta,double deltaV) {
  int i,j;
  double *crddx,*crddy,*crddz;
  double *f;
  struct potential_GOLMAA_dbasin ene,enedx,enedy,enedz;

  crddx=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddy=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddz=(double *)gcemalloc(sizeof(double)*numatom*3);

  f=(double *)gcemalloc(sizeof(double)*3);

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

  GOLMAAff_dbasin_set_calcff(&ene,refcrd1,refcrd2,numatom,R_C_D);
  GOLMAAff_dbasin_set_calcff(&enedx,refcrd1,refcrd2,numatom,R_C_D);
  GOLMAAff_dbasin_set_calcff(&enedy,refcrd1,refcrd2,numatom,R_C_D);
  GOLMAAff_dbasin_set_calcff(&enedz,refcrd1,refcrd2,numatom,R_C_D);

  GOLMAAff_dbasin_calcff(crd,numatom,&ene,delta,deltaV);
  GOLMAAff_dbasin_calcff(crddx,numatom,&enedx,delta,deltaV);
  GOLMAAff_dbasin_calcff(crddy,numatom,&enedy,delta,deltaV);
  GOLMAAff_dbasin_calcff(crddz,numatom,&enedz,delta,deltaV);

  f[0]=-/*0.5**/(enedx.p_t-ene.p_t)/dx*UNIT;
  f[1]=-/*0.5**/(enedy.p_t-ene.p_t)/dx*UNIT;
  f[2]=-/*0.5**/(enedz.p_t-ene.p_t)/dx*UNIT;

  return f;
}
