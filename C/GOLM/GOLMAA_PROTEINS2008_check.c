
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008_check.h"

#include "FFL.h"
#include "TOPO.h"
#include "LA.h"
#include "mymath.h"

#include "PTL.h"

#define UNIT 4.184070*100.0

#define ON 1
#define OFF 0

double GOLMAA_PROTEINS2008_calcff_check(char *inputfilename,char *reffilename,char *parmfilename,
					int numspatom,double dx,
					double ep,int nibnum,double criteria,
					double f_natatt[3],double f_repul[3],
					double f_d[3],double f_a[3],double f_b[3],
					double f_d1[4][3],double f_d2[4][3], int nums,
					double f_as[3][3], int numas) {
  int i,j,k,d;
  double *crd,*crddx,*crddy,*crddz,*refcrd;
  int numatom,numres;
  struct potential e;
  struct force f;
  struct potential_GOLMAA_PROTEINS2008 e_GOLM,e_GOLMdx,e_GOLMdy,e_GOLMdz;

  double x[3];

  double **f1,**f2,p_d,p_dx,p_dy,p_dz;
  double f_d3[4][3];

  double p_natatt,p_natatt_dx,p_natatt_dy,p_natatt_dz;
  double f_natatt3[2][3];
  int atomi,atomj;

  double p_a,p_adx,p_ady,p_adz;
  double f_as0[3][3],f_adx[3][3],f_ady[3][3],f_adz[3][3];

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
  GOLMAA_PROTEINS2008_ff_set_calcff(&e_GOLM,refcrd,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);
  GOLMAA_PROTEINS2008_ff_set_calcff(&e_GOLMdx,refcrd,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);
  GOLMAA_PROTEINS2008_ff_set_calcff(&e_GOLMdy,refcrd,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);
  GOLMAA_PROTEINS2008_ff_set_calcff(&e_GOLMdz,refcrd,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);

  GOLMAA_PROTEINS2008_ff_calcff(crd,numatom,&e_GOLM);
  GOLMAA_PROTEINS2008_ff_calcff(crddx,numatom,&e_GOLMdx);
  GOLMAA_PROTEINS2008_ff_calcff(crddy,numatom,&e_GOLMdy);
  GOLMAA_PROTEINS2008_ff_calcff(crddz,numatom,&e_GOLMdz);

  f_natatt[0]=-(e_GOLMdx.p_natatt_t-e_GOLM.p_natatt_t)/dx*4.184070*100.0;
  f_natatt[1]=-(e_GOLMdy.p_natatt_t-e_GOLM.p_natatt_t)/dx*4.184070*100.0;
  f_natatt[2]=-(e_GOLMdz.p_natatt_t-e_GOLM.p_natatt_t)/dx*4.184070*100.0;

  f_repul[0]=-(e_GOLMdx.p_repul_t-e_GOLM.p_repul_t)/dx*4.184070*100.0;
  f_repul[1]=-(e_GOLMdy.p_repul_t-e_GOLM.p_repul_t)/dx*4.184070*100.0;
  f_repul[2]=-(e_GOLMdz.p_repul_t-e_GOLM.p_repul_t)/dx*4.184070*100.0;

  f_d[0]=-(e_GOLMdx.p_d_t-e_GOLM.p_d_t)/dx*4.184070*100.0;
  f_d[1]=-(e_GOLMdy.p_d_t-e_GOLM.p_d_t)/dx*4.184070*100.0;
  f_d[2]=-(e_GOLMdz.p_d_t-e_GOLM.p_d_t)/dx*4.184070*100.0;

  f_a[0]=-(e_GOLMdx.p_a_t-e_GOLM.p_a_t)/dx*4.184070*100.0;
  f_a[1]=-(e_GOLMdy.p_a_t-e_GOLM.p_a_t)/dx*4.184070*100.0;
  f_a[2]=-(e_GOLMdz.p_a_t-e_GOLM.p_a_t)/dx*4.184070*100.0;

  f_b[0]=-(e_GOLMdx.p_b_t-e_GOLM.p_b_t)/dx*4.184070*100.0;
  f_b[1]=-(e_GOLMdy.p_b_t-e_GOLM.p_b_t)/dx*4.184070*100.0;
  f_b[2]=-(e_GOLMdz.p_b_t-e_GOLM.p_b_t)/dx*4.184070*100.0;

  //////////////////////////////////////////////////////////////////////////////////////////////

  numspatom=(e_GOLM).pairs_angl[numas][0];

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

  GOLMAA_PROTEINS2008_ff_calcANGLE_forcheck(crd,numatom,(e_GOLM).Ka,(e_GOLM).ang_equ,
					    (e_GOLM).pairs_angl,(e_GOLM).num_angl,
					    &p_a,f_as0,numas);
  GOLMAA_PROTEINS2008_ff_calcANGLE_forcheck(crddx,numatom,(e_GOLMdx).Ka,(e_GOLMdx).ang_equ,
					    (e_GOLMdx).pairs_angl,(e_GOLMdx).num_angl,
					    &p_adx,f_adx,numas);
  GOLMAA_PROTEINS2008_ff_calcANGLE_forcheck(crddy,numatom,(e_GOLMdy).Ka,(e_GOLMdy).ang_equ,
					    (e_GOLMdy).pairs_angl,(e_GOLMdy).num_angl,
					    &p_ady,f_ady,numas);
  GOLMAA_PROTEINS2008_ff_calcANGLE_forcheck(crddz,numatom,(e_GOLMdz).Ka,(e_GOLMdz).ang_equ,
					    (e_GOLMdz).pairs_angl,(e_GOLMdz).num_angl,
					    &p_adz,f_adz,numas);

  f_as[0][0]=-(p_adx-p_a)/dx*4.184070*100.0;
  f_as[0][1]=-(p_ady-p_a)/dx*4.184070*100.0;
  f_as[0][2]=-(p_adz-p_a)/dx*4.184070*100.0;

  numspatom=(e_GOLM).pairs_angl[numas][1];

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

  GOLMAA_PROTEINS2008_ff_calcANGLE_forcheck(crd,numatom,(e_GOLM).Ka,(e_GOLM).ang_equ,
					    (e_GOLM).pairs_angl,(e_GOLM).num_angl,
					    &p_a,f_as0,numas);
  GOLMAA_PROTEINS2008_ff_calcANGLE_forcheck(crddx,numatom,(e_GOLMdx).Ka,(e_GOLMdx).ang_equ,
					    (e_GOLMdx).pairs_angl,(e_GOLMdx).num_angl,
					    &p_adx,f_adx,numas);
  GOLMAA_PROTEINS2008_ff_calcANGLE_forcheck(crddy,numatom,(e_GOLMdy).Ka,(e_GOLMdy).ang_equ,
					    (e_GOLMdy).pairs_angl,(e_GOLMdy).num_angl,
					    &p_ady,f_ady,numas);
  GOLMAA_PROTEINS2008_ff_calcANGLE_forcheck(crddz,numatom,(e_GOLMdz).Ka,(e_GOLMdz).ang_equ,
					    (e_GOLMdz).pairs_angl,(e_GOLMdz).num_angl,
					    &p_adz,f_adz,numas);

  f_as[1][0]=-(p_adx-p_a)/dx*4.184070*100.0;
  f_as[1][1]=-(p_ady-p_a)/dx*4.184070*100.0;
  f_as[1][2]=-(p_adz-p_a)/dx*4.184070*100.0;

  numspatom=(e_GOLM).pairs_angl[numas][2];

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

  GOLMAA_PROTEINS2008_ff_calcANGLE_forcheck(crd,numatom,(e_GOLM).Ka,(e_GOLM).ang_equ,
					    (e_GOLM).pairs_angl,(e_GOLM).num_angl,
					    &p_a,f_as0,numas);
  GOLMAA_PROTEINS2008_ff_calcANGLE_forcheck(crddx,numatom,(e_GOLMdx).Ka,(e_GOLMdx).ang_equ,
					    (e_GOLMdx).pairs_angl,(e_GOLMdx).num_angl,
					    &p_adx,f_adx,numas);
  GOLMAA_PROTEINS2008_ff_calcANGLE_forcheck(crddy,numatom,(e_GOLMdy).Ka,(e_GOLMdy).ang_equ,
					    (e_GOLMdy).pairs_angl,(e_GOLMdy).num_angl,
					    &p_ady,f_ady,numas);
  GOLMAA_PROTEINS2008_ff_calcANGLE_forcheck(crddz,numatom,(e_GOLMdz).Ka,(e_GOLMdz).ang_equ,
					    (e_GOLMdz).pairs_angl,(e_GOLMdz).num_angl,
					    &p_adz,f_adz,numas);

  f_as[2][0]=-(p_adx-p_a)/dx*4.184070*100.0;
  f_as[2][1]=-(p_ady-p_a)/dx*4.184070*100.0;
  f_as[2][2]=-(p_adz-p_a)/dx*4.184070*100.0;


  return 0.0;
}

double GOLMAA_PROTEINS2008_ff_calcANGLE_forcheck(double *crd,int numatom,
						 double Ka,double *ang_equ,
						 int **pairs,int numangl,
						 double *p_a,double f_a[3][3],int nums){
  int i,j,k,l;
  int ii,jj,kk;
  double atom[3][3];
  double lenij,lenkj;
  double vij[3],vkj[3];
  double cosijk,angijk;
  double f1,f2;
  double p_a_t=0.0;

  for (i=0;i<3;++i) for (j=0;j<3;++j) f_a[i][j] = 0.0;

  i=nums;
  ii=pairs[i][0];
  jj=pairs[i][1];
  kk=pairs[i][2];
  for (j=0;j<3;++j) {
    atom[0][j]=crd[ii*3+j];
    atom[1][j]=crd[jj*3+j];
    atom[2][j]=crd[kk*3+j];
  }

  lenij = len(atom[0],atom[1]);
  lenkj = len(atom[2],atom[1]);
  for (j=0;j<3;++j) {
    vij[j]=atom[1][j]-atom[0][j];
    vkj[j]=atom[1][j]-atom[2][j];
  }
  cosijk=inprod(vij,vkj,3);
  cosijk=cosijk/lenij/lenkj;
  angijk = acos(cosijk);

  angijk = ang(atom[0],atom[1],atom[2]);

  *p_a = Ka*(angijk-ang_equ[i])*(angijk-ang_equ[i]);
  
  for (j=0;j<3;++j) {
    f1 = -2.0*Ka*(angijk-ang_equ[i])/(lenij*sin(angijk))*(vkj[j]/lenkj-cosijk*vij[j]/lenij)*4.184070*100.0;
    f2 = -2.0*Ka*(angijk-ang_equ[i])/(lenkj*sin(angijk))*(vij[j]/lenij-cosijk*vkj[j]/lenkj)*4.184070*100.0;
    
    f_a[0][j] += f1;
    f_a[1][j] += f2;
    f_a[2][j] += -f1-f2;
  }

  return p_a_t;
}
