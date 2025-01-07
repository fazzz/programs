
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008_check.h"
#include "GOLMAA_MB_PROTEINS2008.h"

#include "FFL.h"
#include "TOPO.h"
#include "LA.h"
#include "mymath.h"

#include "PTL.h"

#define UNIT 4.184070*100.0

#define ON 1
#define OFF 0

double GOLMAA_MB_PROTEINS2008_calcff_check(char *inputfilename,
					   char *reffilename1,char *reffilename2,char *parmfilename,
					   int numspatom,double dx,
					   double d,double de,
					   double ep,int nibnum,double criteria,
					   double f_MB[3],double f_e1[3],double f_e2[3]) {
  int i,j,k;
  double *crd,*crddx,*crddy,*crddz,*refcrd1,*refcrd2;
  int numatom,numres;
  double d2;
  struct potential e;
  struct force f;
  struct potential_GOLMAA_MB_PROTEINS2008 e_GOLM,e_GOLMdx,e_GOLMdy,e_GOLMdz;

  double x[3];

  char *line;
  size_t len=0;

  FILE *inputfile,*reffile1,*reffile2,*outputfile,*parmfile;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  numres=AP.NRES;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd1=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd2=(double *)gcemalloc(sizeof(double)*numatom*3);

  crddx=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddy=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddz=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&i);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  reffile1=efopen(reffilename1,"r");
  getline(&line,&len,reffile1);
  fscanf(reffile1,"%d",&i);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(reffile1,"%lf",&refcrd1[i*3+j]);
  fclose(reffile1);

  reffile2=efopen(reffilename2,"r");
  getline(&line,&len,reffile2);
  fscanf(reffile2,"%d",&i);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(reffile2,"%lf",&refcrd2[i*3+j]);
  fclose(reffile2);

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
  GOLMAA_MB_PROTEINS2008_ff_calcff_set(&e_GOLM,refcrd1,refcrd2,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);
  GOLMAA_MB_PROTEINS2008_ff_calcff_set(&e_GOLMdx,refcrd1,refcrd2,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);
  GOLMAA_MB_PROTEINS2008_ff_calcff_set(&e_GOLMdy,refcrd1,refcrd2,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);
  GOLMAA_MB_PROTEINS2008_ff_calcff_set(&e_GOLMdz,refcrd1,refcrd2,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);

  d2=d*d;
  GOLMAA_MB_PROTEINS2008_ff_calcff(crd,numatom,de,d2,&e_GOLM);
  GOLMAA_MB_PROTEINS2008_ff_calcff(crddx,numatom,de,d2,&e_GOLMdx);
  GOLMAA_MB_PROTEINS2008_ff_calcff(crddy,numatom,de,d2,&e_GOLMdy);
  GOLMAA_MB_PROTEINS2008_ff_calcff(crddz,numatom,de,d2,&e_GOLMdz);

  f_MB[0]=-(e_GOLMdx.p_MB-e_GOLM.p_MB)/dx*4.184070*100.0;
  f_MB[1]=-(e_GOLMdy.p_MB-e_GOLM.p_MB)/dx*4.184070*100.0;
  f_MB[2]=-(e_GOLMdz.p_MB-e_GOLM.p_MB)/dx*4.184070*100.0;

  f_e1[0]=-(e_GOLMdx.e1.p_t-e_GOLM.e1.p_t)/dx*4.184070*100.0;
  f_e1[1]=-(e_GOLMdy.e1.p_t-e_GOLM.e1.p_t)/dx*4.184070*100.0;
  f_e1[2]=-(e_GOLMdz.e1.p_t-e_GOLM.e1.p_t)/dx*4.184070*100.0;

  f_e2[0]=-(e_GOLMdx.e2.p_t-e_GOLM.e2.p_t)/dx*4.184070*100.0;
  f_e2[1]=-(e_GOLMdy.e2.p_t-e_GOLM.e2.p_t)/dx*4.184070*100.0;
  f_e2[2]=-(e_GOLMdz.e2.p_t-e_GOLM.e2.p_t)/dx*4.184070*100.0;

  return 0.0;
}

