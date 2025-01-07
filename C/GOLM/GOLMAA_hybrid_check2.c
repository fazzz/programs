#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLMAA_hybrid_set.h"
#include "GOLMAA_hybrid_check2.h"
#include "FFL.h"

#include "TOPO.h"
#include "LA.h"
#include "mymath.h"

#include "PTL.h"

#define UNIT 4.184070*100.0

#define ON 1
#define OFF 0

double GOLM_AA_hybrid_calcff_check2(char *inputfilename,char *reffilename,char *parmfilename,
				    int numspatom,double dx,
				    double f_natatt[3],double f_repul[3],double f_d[3],double f_a[3],double f_b[3],
				    int nums) {
  int i,j,k,d;
  double *crd,*crddx,*crddy,*crddz,*refcrd;
  int numatom,numres;

  struct potential e,e_dx,e_dy,e_dz;
  struct force f,f_dx,f_dy,f_dz;
  struct potential_GOLMAA_hybrid e_GOLM,e_GOLMdx,e_GOLMdy,e_GOLMdz;

  double x[3];

  double **f1,**f2,p_d,p_dx,p_dy,p_dz;

  char *line;
  size_t len=0;

  FILE *inputfile,*reffile,*outputfile,*parmfile;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  numres=AP.NRES;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);

  crddx=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddy=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddz=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  reffile=efopen(reffilename,"r");
  getline(&line,&len,reffile);
  fscanf(reffile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(reffile,"%lf",&refcrd[i*3+j]);
  fclose(reffile);

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

  GOLMAA_hybrid_ff_set_calcff(&e_GOLM,refcrd,numatom,numres,e.parm.indexnb,e.parm.numnb);
  GOLMAA_hybrid_ff_set_calcff(&e_GOLMdx,refcrd,numatom,numres,e_dx.parm.indexnb,e_dx.parm.numnb);
  GOLMAA_hybrid_ff_set_calcff(&e_GOLMdy,refcrd,numatom,numres,e_dy.parm.indexnb,e_dy.parm.numnb);
  GOLMAA_hybrid_ff_set_calcff(&e_GOLMdz,refcrd,numatom,numres,e_dz.parm.indexnb,e_dz.parm.numnb);

  ffL_calcffandforce_14vdWDAB_woH(crd,numatom,&e,&f);
  ffL_calcffandforce_14vdWDAB_woH(crddx,numatom,&e_dx,&f_dx);
  ffL_calcffandforce_14vdWDAB_woH(crddy,numatom,&e_dy,&f_dy);
  ffL_calcffandforce_14vdWDAB_woH(crddz,numatom,&e_dz,&f_dz);

  GOLMAA_hyb_ff_calcff(crd,numatom,&e_GOLM);
  GOLMAA_hyb_ff_calcff(crddx,numatom,&e_GOLMdx);
  GOLMAA_hyb_ff_calcff(crddy,numatom,&e_GOLMdy);
  GOLMAA_hyb_ff_calcff(crddz,numatom,&e_GOLMdz);

  f_natatt[0]=-(e_GOLMdx.p_natatt_t-e_GOLM.p_natatt_t)/dx*4.184070*100.0;
  f_natatt[1]=-(e_GOLMdy.p_natatt_t-e_GOLM.p_natatt_t)/dx*4.184070*100.0;
  f_natatt[2]=-(e_GOLMdz.p_natatt_t-e_GOLM.p_natatt_t)/dx*4.184070*100.0;

  f_repul[0]=-(e_GOLMdx.p_repul_t-e_GOLM.p_repul_t)/dx*4.184070*100.0;
  f_repul[1]=-(e_GOLMdy.p_repul_t-e_GOLM.p_repul_t)/dx*4.184070*100.0;
  f_repul[2]=-(e_GOLMdz.p_repul_t-e_GOLM.p_repul_t)/dx*4.184070*100.0;

  f_d[0]=-(e_dx.p_d_t-e.p_d_t)/dx*4.184070*100.0;
  f_d[1]=-(e_dy.p_d_t-e.p_d_t)/dx*4.184070*100.0;
  f_d[2]=-(e_dz.p_d_t-e.p_d_t)/dx*4.184070*100.0;

  f_a[0]=-(e_dx.p_a_t-e.p_a_t)/dx*4.184070*100.0;
  f_a[1]=-(e_dy.p_a_t-e.p_a_t)/dx*4.184070*100.0;
  f_a[2]=-(e_dz.p_a_t-e.p_a_t)/dx*4.184070*100.0;

  f_b[0]=-(e_dx.p_b_t-e.p_b_t)/dx*4.184070*100.0;
  f_b[1]=-(e_dy.p_b_t-e.p_b_t)/dx*4.184070*100.0;
  f_b[2]=-(e_dz.p_b_t-e.p_b_t)/dx*4.184070*100.0;

  return 0;

}
