
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "FFL.h"
#include "TOPO.h"
#include "LA.h"
#include "mymath.h"

#include "PTL.h"

#define UNIT 4.184070*100.0

#define ON 1
#define OFF 0

double AAFF_Amber_calcff_check(char *inputfilename,char *parmfilename,
			       int numspatom,double dx,
			       double f_e[3],double f_LJ[3],
			       double f_e_14[3],double f_LJ_14[3],
			       double f_d[3],double f_a[3],double f_b[3]) {
  int i,j,k,d;
  double *crd,*crddx,*crddy,*crddz;
  int numatom,numres;
  struct potential e,edx,edy,edz;
  struct force f,fdx,fdy,fdz;

  char *line;
  size_t len=0;

  FILE *inputfile,*outputfile,*parmfile;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  numres=AP.NRES;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  crddx=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddy=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddz=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

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
  ffL_set_calcffandforce(&edx,&fdx);
  ffL_set_calcffandforce(&edy,&fdy);
  ffL_set_calcffandforce(&edz,&fdz);

  ffL_calcffandforce(crd,numatom,&e,&f);
  ffL_calcffandforce(crddx,numatom,&edx,&fdx);
  ffL_calcffandforce(crddy,numatom,&edy,&fdy);
  ffL_calcffandforce(crddz,numatom,&edz,&fdz);

  f_e[0]=-(edx.p_e_t-e.p_e_t)/dx*4.184070*100.0;
  f_e[1]=-(edy.p_e_t-e.p_e_t)/dx*4.184070*100.0;
  f_e[2]=-(edz.p_e_t-e.p_e_t)/dx*4.184070*100.0;

  f_LJ[0]=-(edx.p_LJ_t-e.p_LJ_t)/dx*4.184070*100.0;
  f_LJ[1]=-(edy.p_LJ_t-e.p_LJ_t)/dx*4.184070*100.0;
  f_LJ[2]=-(edz.p_LJ_t-e.p_LJ_t)/dx*4.184070*100.0;

  f_e_14[0]=-(edx.p_e_14_t-e.p_e_14_t)/dx*4.184070*100.0;
  f_e_14[1]=-(edy.p_e_14_t-e.p_e_14_t)/dx*4.184070*100.0;
  f_e_14[2]=-(edz.p_e_14_t-e.p_e_14_t)/dx*4.184070*100.0;

  f_LJ_14[0]=-(edx.p_LJ_14_t-e.p_LJ_14_t)/dx*4.184070*100.0;
  f_LJ_14[1]=-(edy.p_LJ_14_t-e.p_LJ_14_t)/dx*4.184070*100.0;
  f_LJ_14[2]=-(edz.p_LJ_14_t-e.p_LJ_14_t)/dx*4.184070*100.0;

  f_d[0]=-(edx.p_d_t-e.p_d_t)/dx*4.184070*100.0;
  f_d[1]=-(edy.p_d_t-e.p_d_t)/dx*4.184070*100.0;
  f_d[2]=-(edz.p_d_t-e.p_d_t)/dx*4.184070*100.0;

  f_a[0]=-(edx.p_a_t-e.p_a_t)/dx*4.184070*100.0;
  f_a[1]=-(edy.p_a_t-e.p_a_t)/dx*4.184070*100.0;
  f_a[2]=-(edz.p_a_t-e.p_a_t)/dx*4.184070*100.0;

  f_b[0]=-(edx.p_b_t-e.p_b_t)/dx*4.184070*100.0;
  f_b[1]=-(edy.p_b_t-e.p_b_t)/dx*4.184070*100.0;
  f_b[2]=-(edz.p_b_t-e.p_b_t)/dx*4.184070*100.0;

  return 0.0;
}

