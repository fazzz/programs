#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLMAA_hybrid_set.h"
#include "GOLMAA_hybrid_check.h"
#include "FFL.h"
#include "PTL.h"

#define ON 1
#define OFF 0

double MD_GOLMAAff_hybrid_calcff_check(double *crd, int numatom,
				       FILE *parmtop;
				       int numspatom,double dx,
				       double f_b[3],double f_a[3],double f_d[3],
				       double f_e_14[3],double f_LJ_14[3]) {
  int i,j;
  double *crddx,*crddy,*crddz;
  struct potential e,e_dx,e_dy,e_dz;
  struct force f,f_dx,f_dy,f_dz;

  readParmtopL(parmfile);

  crddx=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddy=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddz=(double *)gcemalloc(sizeof(double)*numatom*3);

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

  ffL_set_calcffandforce(&e,&f);
  ffL_set_calcffandforce(&e_dx,&f_dx);
  ffL_set_calcffandforce(&e_dy,&f_dy);
  ffL_set_calcffandforce(&e_dz,&f_dz);

  ffL_calcffandforce(crd,numatom,&e,&f);
  ffL_calcffandforce(crddx,numatom,&e_dx,&f_dx);
  ffL_calcffandforce(crddy,numatom,&e_dy,&f_dy);
  ffL_calcffandforce(crddz,numatom,&e_dz,&f_dz);

  f_b[0]=-(e_dx.p_b[numspatom]-e.p_b[numspatom])/dx*UNIT;
  f_b[1]=-(e_dy.p_b[numspatom]-e.p_b[numspatom])/dx*UNIT;
  f_b[2]=-(e_dz.p_b[numspatom]-e.p_b[numspatom])/dx*UNIT;

  f_a[0]=-(e_dx.p_a[numspatom]-e.p_a[numspatom])/dx*UNIT;
  f_a[1]=-(e_dy.p_a[numspatom]-e.p_a[numspatom])/dx*UNIT;
  f_a[2]=-(e_dz.p_a[numspatom]-e.p_a[numspatom])/dx*UNIT;

  f_d[0]=-(e_dx.p_d[numspatom]-e.p_d[numspatom])/dx*UNIT;
  f_d[1]=-(e_dy.p_d[numspatom]-e.p_d[numspatom])/dx*UNIT;
  f_d[2]=-(e_dz.p_d[numspatom]-e.p_d[numspatom])/dx*UNIT;

  f_e_14[0]=-(e_dx.p_e_14[numspatom]-e.p_e_14[numspatom])/dx*UNIT;
  f_e_14[1]=-(e_dy.p_e_14[numspatom]-e.p_e_14[numspatom])/dx*UNIT;
  f_e_14[2]=-(e_dz.p_e_14[numspatom]-e.p_e_14[numspatom])/dx*UNIT;

  f_LJ_14[0]=-(e_dx.p_LJ_14[numspatom]-e.p_LJ_14[numspatom])/dx*UNIT;
  f_LJ_14[1]=-(e_dy.p_LJ_14[numspatom]-e.p_LJ_14[numspatom])/dx*UNIT;
  f_LJ_14[2]=-(e_dz.p_LJ_14[numspatom]-e.p_LJ_14[numspatom])/dx*UNIT;

}
