
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "FF.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

int main(int argc, char *argv[]) {
  int i,j,k,flag;

  double *ele,*ALJ,*BLJ;
  double *p_e,*p_LJ,p_t,p_et,p_LJt;
  double *f_e,*f_LJ,f_t;
  int numnb, *indexnb;
  int numatom,numpara;
  double *crd;
  char *line;
  size_t len=0;

  char *inputfilename,*inputfilename2;
  char *outputfilename;
  FILE *inputfile,*inputfile2;
  FILE *outputfile;

  if (argc < 3) {
    printf("USAGE: ./CFF.exe inputfilename(crd) inputfilename2(top) \n");
    exit(1);
  }
  inputfilename   = *++argv;
  inputfilename2  = *++argv;
  outputfilename  = *++argv;

  inputfile2=efopen(inputfilename2,"r");
  readParmtop(inputfile2);
  fclose(inputfile2);
  numatom=AP.NATOM;
  numpara=AP.NTYPES*(AP.NTYPES+1)/2;
  numnb=ff_set_numnb();

  indexnb=(int *)gcemalloc(sizeof(int)*numnb);
  ele=(double *)gcemalloc(sizeof(double)*numatom);
  ALJ=(double *)gcemalloc(sizeof(double)*numpara);
  BLJ=(double *)gcemalloc(sizeof(double)*numpara);
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  p_e=(double *)gcemalloc(sizeof(double)*numnb);
  p_LJ=(double *)gcemalloc(sizeof(double)*numnb);

  ff_set_NB_PARM(ele,ALJ,BLJ,numatom);
  ff_set_NB_index(indexnb,numnb,numatom);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  getline(&line,&len,inputfile);
  io_scanconf(inputfile,numatom,crd,'x');
  fclose(inputfile);

  ff_calcFFNB(ele,ALJ,BLJ,p_e,p_LJ,f_e,f_LJ,numnb,indexnb,numatom,crd,2,0);

  p_et=0.0;p_LJt=0.0;
  for (i=0;i<numnb;++i) {
    p_et+=p_e[i];
    p_LJt+=p_LJ[i];
  }
  p_t=p_et+p_LJt;
  
  printf("/***********************************************/\n");
  printf("poe_energy  = %e kcal/mol  \n",p_t);
  printf("ele_energy  = %e kcal/mol  \n",p_et);
  printf("VDW_energy  = %e kcal/mol  \n",p_LJt);
  printf("/***********************************************/\n");

  return 0;
}

