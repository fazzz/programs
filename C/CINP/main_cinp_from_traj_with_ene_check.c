
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "CINP.h"
#include "FF.h"
#include "PT.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j,amberflag;
  int numsnap,numinpcrd,numinterval;
  int *num;
  int f;
  double checkv;
  double *crd;

  double *ele,*ALJ,*BLJ;
  double *p_e,*p_LJ,*p_d,*p_t;
  double *f_e,*f_LJ,*n_d;
  double *p_e_14,*p_LJ_14;
  double *f_e_14,*f_LJ_14;
  int numnb,num14, *indexnb,*index14;
  int numatom,numpara;

  char *option,*opt,*dummy;  
  char *inputfilename1,*inputfilename2,*inputfilename3,*outputfilenamebase,*ext;
  FILE *inputfile1,*inputfile2,*inputfile3;
  
  if (argc < 6) {
    printf("USAGE: ./%s -[at] numsnap numinterval checkv inputfilename1(trj) inputfilename2(parm) inputfilename3(index)  outputfilenamebase(inp) \n",argv[0]);
    exit(1);
  }
  amberflag=(*++argv)[0];
  if (amberflag=='a')
    amberflag=ON;
  else
    amberflag=OFF;
 
  numsnap=atoi(*++argv);
  numinterval=atoi(*++argv);
  checkv=atof(*++argv);
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  inputfilename3 = *++argv;
  outputfilenamebase = *++argv;
  
  numinpcrd=(int)numsnap/numinterval;
  num=(int *)gcemalloc(sizeof(int)*numinpcrd);
  for (i=0;i<numinpcrd;++i)
    num[i]=i*numinterval;
  
  inputfile2=efopen(inputfilename2,"r");
  readParmtop(inputfile2);
  fclose(inputfile2);
  numatom=AP.NATOM;

  indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  ele=(double *)gcemalloc(sizeof(double)*numatom);
  ALJ=(double *)gcemalloc(sizeof(double)*numpara);
  BLJ=(double *)gcemalloc(sizeof(double)*numpara);
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  p_e=(double *)gcemalloc(sizeof(double)*numnb);
  p_LJ=(double *)gcemalloc(sizeof(double)*numnb);
  index14=(int *)gcemalloc(sizeof(int)*num14*2);
  p_e_14=(double *)gcemalloc(sizeof(double)*numnb);
  p_LJ_14=(double *)gcemalloc(sizeof(double)*numnb);
  p_d=(double *)gcemalloc(sizeof(double)*(AP.NPHIH+AP.MPHIA));

  ff_set_NB_PARM(ele,ALJ,BLJ,numatom);

  inputfile3=efopen(inputfilename3,"r");
  for (i=0;i<numnb;++i) {
    fscanf(inputfile3,"%d",&f);indexnb[i*2]=f-1;
    fscanf(inputfile3,"%s",&dummy);
    fscanf(inputfile3,"%d",&f);indexnb[i*2+1]=f-1;
  }
  for (i=0;i<num14;++i) {
    fscanf(inputfile3,"%d",&f);index14[i*2]=f-1;
    fscanf(inputfile3,"%s",&dummy);
    fscanf(inputfile3,"%d",&f);index14[i*2+1]=f-1;
  }
  fclose(inputfile3);

  inputfile1=efopen(inputfilename1,"r");

  p_t=(double *)gcemalloc(sizeof(double)*numsnap);
  for (i=0;i<numsnap;++i) {
    io_scanconf(inputfile1,numatom,crd,'x');
    ff_calcFFNB(ele,ALJ,BLJ,p_e,p_LJ,f_e,f_LJ,numnb,indexnb,numatom,crd,2,0);
    ff_calcFFNB(ele,ALJ,BLJ,p_e_14,p_LJ_14,f_e_14,f_LJ_14,num14,index14,numatom,crd,2,0);
    ff_calcDIHE(p_d,n_d,crd,1,0,0);

    p_t[i]=0.0;
    for (j=0;j<numnb;++j)
      p_t[i]+=p_e[j]+p_LJ[j];
    for (j=0;j<num14;++j)
      p_t[i]+=1.0/1.2*p_e_14[j]+0.5*p_LJ_14[j];
    for (j=0;j<AP.NPHIH+AP.MPHIA;++j)
      p_t[i]+=p_d[j];
  }
  fclose(inputfile1);

  //  inputfile1=efopen(inputfilename1,"r");  
  Create_inpcrd_from_trj_with_ene_check(outputfilenamebase,inputfilename1,numsnap,numinpcrd,numatom,num,amberflag,p_t,checkv);
  //  fclose(inputfile1);  
  
  return 0;
}


