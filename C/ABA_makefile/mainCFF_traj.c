
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "FF.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0

int main(int argc, char *argv[]) {
  int i,j,k,flag,flagt,amberflag;
  int f;

  double *ele,*ALJ,*BLJ;
  double *p_e,*p_LJ,*p_d,*p_a,*p_b,p_t;
  double p_e_t,p_LJ_t,p_e_14_t,p_LJ_14_t,p_d_t,p_a_t,p_b_t;
  double *f_e,*f_LJ,*n_d;
  double *p_e_14,*p_LJ_14;
  double *f_e_14,*f_LJ_14;
  int numnb,num14, *indexnb,*index14;
  int numatom,numpara,numstep;
  double *crd;
  char *line,*dummy;
  size_t len=0;

  char *inputfilename,*inputfilename2,*inputfilename3;
  char *outputfilename,*outputfilename2,*outputfilename3,*outputfilename4,*outputfilename5,*outputfilename6,*outputfilenamea,*outputfilenameb;
  FILE *inputfile,*inputfile2,*inputfile3;
  FILE *outputfile,*outputfile2,*outputfile3,*outputfile4,*outputfile5,*outputfile6,*outputfilea,*outputfileb;

  if (argc < 5) {
    printf("USAGE: ./%s [at]  [54db] [to] numnb num14 numstep inputfilename(crd) inputfilename2(top) inputfilename3(index) outputfilename(ene_es) outputfilename2(ene_LJ) outputfilename3(ene_es_1-4) outputfilename4(ene_LJ_1-4) outputfilename5(ene_di) outputfilenamea(ene_ang) outputfilenameb(ene_bon) outputfilename6(ene_t)\n",argv[0]);
    exit(1);
  }
  flag=(*++argv)[0];
  if (flag != 'a' && flag != 't') {
    printf("flag error: must be a or 4 or t");
    exit(1);
  }  
  if (flag=='a')
    amberflag=ON;
  else
    amberflag=OFF;
  flag=(*++argv)[0];
  if (flag != '5' && flag != '4' && flag != 'd' && flag != 'b') {
    printf("flag error: must be 5 or 4 or d or b");
    exit(1);
  }
  flagt=(*++argv)[0];
  if (flagt != 't' && flagt != 'o') {
    printf("flag error: must be t or o ");
    exit(1);
  }
  numnb=atoi(*++argv);
  num14=atoi(*++argv);
  numstep=atoi(*++argv);
  inputfilename   = *++argv;
  inputfilename2  = *++argv;
  inputfilename3  = *++argv;
  outputfilename  = *++argv;
  outputfilename2  = *++argv;
  outputfilename3  = *++argv;
  outputfilename4  = *++argv;
  outputfilename5  = *++argv;
  outputfilenamea  = *++argv;
  outputfilenameb  = *++argv;
  outputfilename6  = *++argv;

  inputfile2=efopen(inputfilename2,"r");
  readParmtop(inputfile2);
  fclose(inputfile2);
  numatom=AP.NATOM;
  numpara=AP.NTYPES*(AP.NTYPES+1)/2;

  indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  ele=(double *)gcemalloc(sizeof(double)*numatom);
  ALJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);
  BLJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  p_e=(double *)gcemalloc(sizeof(double)*numnb);
  p_LJ=(double *)gcemalloc(sizeof(double)*numnb);
  index14=(int *)gcemalloc(sizeof(int)*num14*2);
  p_e_14=(double *)gcemalloc(sizeof(double)*numnb);
  p_LJ_14=(double *)gcemalloc(sizeof(double)*numnb);
  p_d=(double *)gcemalloc(sizeof(double)*(AP.NPHIH+AP.MPHIA));
  p_a=(double *)gcemalloc(sizeof(double)*(AP.NTHETH+AP.MTHETA));
  p_b=(double *)gcemalloc(sizeof(double)*(AP.NBONH+AP.MBONA));

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

  inputfile=efopen(inputfilename,"r");
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");
  if (flag=='4' || flag=='d' || flag=='b' ) {
    outputfile3=efopen(outputfilename3,"w");
    outputfile4=efopen(outputfilename4,"w");
  }
  if (flag=='d' || flag=='b' ) {
    outputfile5=efopen(outputfilename5,"w");
  }
  if (flag=='b' ) {
    outputfilea=efopen(outputfilenamea,"w");
    outputfileb=efopen(outputfilenameb,"w");
  }
  outputfile6=efopen(outputfilename6,"w");

  if (amberflag==ON)
    getline(&line,&len,inputfile);
  for (i=0;i<numstep;++i) {
    io_scanconf(inputfile,numatom,crd,'x');
    ff_calcFFNB(ele,ALJ,BLJ,p_e,p_LJ,f_e,f_LJ,numnb,indexnb,numatom,crd,2,0);
    if (flag=='4' || flag=='d' || flag=='b')
      ff_calcFFNB(ele,ALJ,BLJ,p_e_14,p_LJ_14,f_e_14,f_LJ_14,num14,index14,numatom,crd,2,0);
    if (flag=='d' || flag=='b')
      ff_calcDIHE(p_d,n_d,crd,1,0,0);
    if (flag=='b') {
      ff_calcANGLE(p_a,crd);
      ff_calcBOND(p_b,crd);
    }

    p_t=0.0;
    p_e_t=0.0;
    p_LJ_t=0.0;
    p_e_14_t=0.0;
    p_LJ_14_t=0.0;
    p_d_t=0.0;
    p_a_t=0.0;
    p_b_t=0.0;
    
    for (j=0;j<numnb;++j) {
      p_t+=p_e[j]+p_LJ[j];
      if (flagt=='t') {
	p_e_t+=p_e[j];
	p_LJ_t+=p_LJ[j];
      }
    }
    if (flag=='4' || flag=='d' || flag=='b')
      for (j=0;j<num14;++j) {
	p_t+=1.0/1.2*p_e_14[j]+0.5*p_LJ_14[j];
	if (flagt=='t') {
	  p_e_14_t+=p_e_14[j];
	  p_LJ_14_t+=p_LJ_14[j];
	}
      }
    if (flag=='d' || flag=='b') {
      for (j=0;j<AP.NPHIH+AP.MPHIA;++j) {
	p_t+=p_d[j];
	if (flagt=='t') {
	  p_d_t+=p_d[j];
	}
      }
    }
    if (flag=='b') {
      for (j=0;j<AP.NTHETH+AP.MTHETA;++j) {
	p_t+=p_a[j];
	p_t+=p_b[j];
	if (flagt=='t') {
	  p_a_t+=p_a[j];
	  p_b_t+=p_b[j];
	}
      }
    }

    fprintf(outputfile,"%d ",i);
    if (flagt=='t')
      fprintf(outputfile,"%e \n",p_e_t);
    else {
      for (j=0;j<numnb;++j)
	fprintf(outputfile,"%e ",p_e[j]);
      fprintf(outputfile,"\n");
    }
    fprintf(outputfile2,"%d ",i);
    if (flagt=='t')
      fprintf(outputfile2,"%e \n",p_LJ_t);
    else {
      for (j=0;j<numnb;++j)
	fprintf(outputfile2,"%e ",p_LJ[j]);
      fprintf(outputfile2,"\n");
    }
    if (flag=='4' || flag=='d' || flag=='b') {
      fprintf(outputfile3,"%d ",i);
      if (flagt=='t')
	fprintf(outputfile3,"%e \n",p_e_14_t);
      else {
	for (j=0;j<num14;++j)
	  fprintf(outputfile3,"%e ",1.0/1.2*p_e_14[j]);
	fprintf(outputfile3,"\n");
      }
      fprintf(outputfile4,"%d ",i);
      if (flagt=='t')
	fprintf(outputfile4,"%e \n",p_LJ_14_t);
      else {
	for (j=0;j<num14;++j)
	  fprintf(outputfile4,"%e \n",0.5*p_LJ_14[j]);
	fprintf(outputfile4,"\n");
	fprintf(outputfile5,"%d ",i);
      }
    }
    if (flag=='d' || flag=='b') {
      fprintf(outputfile5,"%d ",i);
      if (flagt=='t')
	fprintf(outputfile5,"%e \n",p_d_t);
      else {
	for (j=0;j<AP.NPHIH+AP.MPHIA;++j)
	  fprintf(outputfile5,"%e ",p_d[j]);
	fprintf(outputfile5,"\n");
      }
    }
    if (flag=='b') {
      fprintf(outputfilea,"%d ",i);
      if (flagt=='t')
	fprintf(outputfilea,"%e \n",p_a_t);
      else {
	for (j=0;j<AP.NTHETH+AP.MTHETA;++j)
	  fprintf(outputfilea,"%e ",p_a[j]);
	fprintf(outputfilea,"\n");
      }
      if (flagt=='t')
	fprintf(outputfileb,"%e \n",p_b_t);
      else {
	for (j=0;j<AP.NBONH+AP.MBONA;++j)
	  fprintf(outputfileb,"%e ",p_b[j]);
	fprintf(outputfileb,"\n");
      }
    }
    fprintf(outputfile6,"%e\n",p_t);
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
  if (flag=='b') {
    fclose(outputfilea);
    fclose(outputfileb);
  }

  fclose(outputfile6);

  return 0;
}

