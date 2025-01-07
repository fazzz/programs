#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLMAA_hybrid_set.h"
#include "GOLMAA_hybrid_check.h"

#define ON 1
#define OFF 0

double GOLMAAff_hybrid_calcff_check(double *crd, int numatom,
				    struct potential_GOLMAA_hybrid *ene,
				    int numspatom,double dx,
				    double f_natatt[3],double f_repul[3]) {
  int i,j;
  double *crddx,*crddy,*crddz;
  double *p_natatt,*p_repul;
  double *p_natatt_dx,*p_natatt_dy,*p_natatt_dz;
  double *p_repul_dx,*p_repul_dy,*p_repul_dz;
  double **f1,**f2;

  crddx=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddy=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddz=(double *)gcemalloc(sizeof(double)*numatom*3);

  p_natatt=(double *)gcemalloc(sizeof(double)*numatom);
  p_natatt_dx=(double *)gcemalloc(sizeof(double)*numatom);
  p_natatt_dy=(double *)gcemalloc(sizeof(double)*numatom);
  p_natatt_dz=(double *)gcemalloc(sizeof(double)*numatom);

  p_repul=(double *)gcemalloc(sizeof(double)*numatom);
  p_repul_dx=(double *)gcemalloc(sizeof(double)*numatom);
  p_repul_dy=(double *)gcemalloc(sizeof(double)*numatom);
  p_repul_dz=(double *)gcemalloc(sizeof(double)*numatom);

  f1=(double **)gcemalloc(sizeof(double *)*numatom);
  f2=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    f1[i]=(double *)gcemalloc(sizeof(double)*3);
    f2[i]=(double *)gcemalloc(sizeof(double)*3);
  }

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

  GOLMAA_hyb_ff_calcff_nonlocal(crd,numatom,
				(*ene).NC_index,(*ene).numNC,(*ene).NotNC_index,(*ene).numNotNC,
				(*ene).ALJ_natatt,(*ene).BLJ_natatt,(*ene).ep_natatt,(*ene).ALJ_repul,
				p_natatt,p_repul,f1,f2);

  GOLMAA_hyb_ff_calcff_nonlocal(crddx,numatom,
				(*ene).NC_index,(*ene).numNC,(*ene).NotNC_index,(*ene).numNotNC,
				(*ene).ALJ_natatt,(*ene).BLJ_natatt,(*ene).ep_natatt,(*ene).ALJ_repul,
				p_natatt_dx,p_repul_dx,f1,f2);

  GOLMAA_hyb_ff_calcff_nonlocal(crddy,numatom,
				(*ene).NC_index,(*ene).numNC,(*ene).NotNC_index,(*ene).numNotNC,
				(*ene).ALJ_natatt,(*ene).BLJ_natatt,(*ene).ep_natatt,(*ene).ALJ_repul,
				p_natatt_dy,p_repul_dy,f1,f2);

  GOLMAA_hyb_ff_calcff_nonlocal(crddz,numatom,
				(*ene).NC_index,(*ene).numNC,(*ene).NotNC_index,(*ene).numNotNC,
				(*ene).ALJ_natatt,(*ene).BLJ_natatt,(*ene).ep_natatt,(*ene).ALJ_repul,
				p_natatt_dz,p_repul_dz,f1,f2);

  f_natatt[0]=-(p_natatt_dx[numspatom]-p_natatt[numspatom])/dx*UNIT;
  f_natatt[1]=-(p_natatt_dy[numspatom]-p_natatt[numspatom])/dx*UNIT;
  f_natatt[2]=-(p_natatt_dz[numspatom]-p_natatt[numspatom])/dx*UNIT;

  f_repul[0]=-(p_repul_dx[numspatom]-p_repul[numspatom])/dx*UNIT;
  f_repul[1]=-(p_repul_dy[numspatom]-p_repul[numspatom])/dx*UNIT;
  f_repul[2]=-(p_repul_dz[numspatom]-p_repul[numspatom])/dx*UNIT;

}
