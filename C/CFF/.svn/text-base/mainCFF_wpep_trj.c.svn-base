
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
  int f;

  double *ele,*ALJ,*BLJ;
  double *p_e,*p_LJ,*p_d,p_t;
  double *f_e,*f_LJ,*n_d;
  double *p_e_14,*p_LJ_14;
  double *f_e_14,*f_LJ_14;
  int numnb,num14, *indexnb,*index14;
  int numatom,numpara,numstep;
  double *u_1_5,*u_1_4,*u_d;
  double fact,dummy2;
  double *crd;
  char *line,*dummy;
  size_t len=0;

  char *inputfilename,*inputfilename2,*inputfilename3,*inputfilename4;
  char *outputfilename,*outputfilename2,*outputfilename3,*outputfilename4,*outputfilename5,*outputfilename6;
  FILE *inputfile,*inputfile2,*inputfile3,*inputfile4;
  FILE *outputfile,*outputfile2,*outputfile3,*outputfile4,*outputfile5,*outputfile6;

  if (argc < 13) {
    printf("USAGE: ./%s [54d] numnb num14 numstep fact inputfilename(crd) inputfilename2(top) inputfilename3(index) inputfilename4(eig) outputfilename(ene_es) outputfilename2(ene_LJ) outputfilename(ene_es_1-4) outputfilename2(ene_LJ_1-4) outputfilename5(ene_t)\n",argv[0]);
    exit(1);
  }
  flag=(*++argv)[0];
  if (flag != '5' && flag != '4' && flag != 'd') {
    printf("flag error: must be 5 or 4 or d");
    exit(1);
  }
  numnb=atoi(*++argv);
  num14=atoi(*++argv);
  numstep=atoi(*++argv);
  fact=atof(*++argv);
  inputfilename   = *++argv;
  inputfilename2  = *++argv;
  inputfilename3  = *++argv;
  inputfilename4  = *++argv;
  outputfilename  = *++argv;
  outputfilename2  = *++argv;
  if (flag=='4' || flag=='d') {
    outputfilename3  = *++argv;
    outputfilename4  = *++argv;
  }
  if (flag=='d') {
    outputfilename5  = *++argv;
  }
  outputfilename6  = *++argv;

  inputfile2=efopen(inputfilename2,"r");
  readParmtop(inputfile2);
  fclose(inputfile2);
  numatom=AP.NATOM;
  numpara=AP.NTYPES*(AP.NTYPES+1)/2;

  u_1_5=(double *)gcemalloc(sizeof(double)*numnb*2);
  if (flag=='4' || flag=='d')
    u_1_4=(double *)gcemalloc(sizeof(double)*num14*2);
  if (flag=='d')
    u_d=(double *)gcemalloc(sizeof(double)*(AP.NPHIH+AP.MPHIA));

  inputfile4=efopen(inputfilename4,"r");
  for (i=0;i<numnb;++i)
    fscanf(inputfile4,"%lf",&u_1_5[i*2]);
  for (i=0;i<numnb;++i)
    fscanf(inputfile4,"%lf",&u_1_5[i*2+1]);
  if (flag=='4' || flag=='d') {
    for (i=0;i<num14;++i)
      fscanf(inputfile4,"%lf",&u_1_4[i*2]);
    for (i=0;i<num14;++i)
      fscanf(inputfile4,"%lf",&u_1_4[i*2+1]);
  }
  if (flag=='d') {
    for (i=0;i<AP.NPHIH;++i)
      fscanf(inputfile4,"%lf",&u_d[i]);
    for (i=0;i<AP.MPHIA;++i)
      fscanf(inputfile4,"%lf",&u_d[i+AP.NPHIH]);
  }
  fclose(inputfile4);

  for (i=0;i<numnb;++i)
    printf("%lf\n",u_1_5[i*2]);
  for (i=0;i<numnb;++i)
    printf("%lf\n",u_1_5[i*2+1]);
  if (flag=='4' || flag=='d') {
    for (i=0;i<num14;++i)
      printf("%lf\n",u_1_4[i*2]);
    for (i=0;i<num14;++i)
      printf("%lf\n",u_1_4[i*2+1]);
  }
  if (flag=='d') {
    for (i=0;i<AP.NPHIH;++i)
      printf("%lf\n",u_d[i]);
    for (i=0;i<AP.MPHIA;++i)
      printf("%lf\n",u_d[i+AP.NPHIH]);
  }
  printf("%lf\n",fact);

  indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  ele=(double *)gcemalloc(sizeof(double)*numatom);
  ALJ=(double *)gcemalloc(sizeof(double)*numpara);
  BLJ=(double *)gcemalloc(sizeof(double)*numpara);
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  p_e=(double *)gcemalloc(sizeof(double)*numnb);
  p_LJ=(double *)gcemalloc(sizeof(double)*numnb);
  if (flag=='4' || flag=='d') {
    index14=(int *)gcemalloc(sizeof(int)*num14*2);
    p_e_14=(double *)gcemalloc(sizeof(double)*numnb);
    p_LJ_14=(double *)gcemalloc(sizeof(double)*numnb);
  }
  if (flag=='d') {
    p_d=(double *)gcemalloc(sizeof(double)*(AP.NPHIH+AP.MPHIA));
  }

  ff_set_NB_PARM(ele,ALJ,BLJ,numatom);

  inputfile3=efopen(inputfilename3,"r");
  for (i=0;i<numnb;++i) {
    fscanf(inputfile3,"%d",&f);indexnb[i*2]=f-1;
    fscanf(inputfile3,"%s",&dummy);
    fscanf(inputfile3,"%d",&f);indexnb[i*2+1]=f-1;
  }
  if (flag=='4' || flag=='d') {
    for (i=0;i<num14;++i) {
      fscanf(inputfile3,"%d",&f);index14[i*2]=f-1;
      fscanf(inputfile3,"%s",&dummy);
      fscanf(inputfile3,"%d",&f);index14[i*2+1]=f-1;
    }
  }
  fclose(inputfile3);

  inputfile=efopen(inputfilename,"r");
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");
  if (flag=='4' || flag=='d') {
    outputfile3=efopen(outputfilename3,"w");
    outputfile4=efopen(outputfilename4,"w");
  }
  if (flag=='d') {
    outputfile5=efopen(outputfilename5,"w");
  }
  outputfile6=efopen(outputfilename6,"w");

  for (i=0;i<numstep;++i) {
    io_scanconf(inputfile,numatom,crd,'x');    
    ff_calcFFNB_wpep(ele,ALJ,BLJ,p_e,p_LJ,f_e,f_LJ,numnb,indexnb,numatom,crd,u_1_5,fact,2,0);
    if (flag=='4' || flag=='d') {
      ff_calcFFNB_wpep(ele,ALJ,BLJ,p_e_14,p_LJ_14,f_e_14,f_LJ_14,num14,index14,numatom,crd,u_1_4,fact,2,0);
    }
    if (flag=='d')
      ff_calcDIHE_wpep(p_d,n_d,crd,u_d,1,0,0);
    p_t=0.0;
    for (j=0;j<numnb;++j)
      p_t+=p_e[j]+p_LJ[j];
    if (flag=='4' || flag=='d')
      for (j=0;j<num14;++j)
	p_t+=1.0/1.2*p_e_14[j]+0.5*p_LJ_14[j];
    if (flag=='d')
      for (j=0;j<AP.NPHIH+AP.MPHIA;++j)
	p_t+=p_d[j];

    fprintf(outputfile,"%d ",i);
    for (j=0;j<numnb;++j)
      fprintf(outputfile,"%e ",p_e[j]);
    fprintf(outputfile,"\n");
    fprintf(outputfile2,"%d ",i);
    for (j=0;j<numnb;++j)
      fprintf(outputfile2,"%e ",p_LJ[j]);
    fprintf(outputfile2,"\n");
    fprintf(outputfile3,"%d ",i);
    if (flag=='4' || flag=='d') {
      for (j=0;j<num14;++j)
	fprintf(outputfile3,"%e ",1.0/1.2*p_e_14[j]);
      fprintf(outputfile3,"\n");
      fprintf(outputfile4,"%d ",i);
      for (j=0;j<num14;++j)
	fprintf(outputfile4,"%e ",0.5*p_LJ_14[j]);
      fprintf(outputfile4,"\n");
    }
    if (flag=='d') {
      for (j=0;j<AP.NPHIH+AP.MPHIA;++j)
	fprintf(outputfile5,"%e ",p_d[j]);
      fprintf(outputfile5,"\n");
    }
    fprintf(outputfile6,"%e \n",p_t);
  }

  fclose(inputfile);
  fclose(outputfile);
  fclose(outputfile2);
  if (flag=='4' || flag=='d') {
    fclose(outputfile3);
    fclose(outputfile4);
  }
  if (flag=='d') {
    fclose(outputfile5);
  }
  fclose(outputfile6);

  return 0;
}

