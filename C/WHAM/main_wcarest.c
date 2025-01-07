
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "WHAM.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

int main(int argc, char *argv[]) {
  int i,j,k;
  int AMBERFLAG=OFF,flag;
  char *line;
  size_t len=0;

  int nwindows,numatom,numca,*nt;
  int nttemp;
  double *ebF,***ebW;
  double **crd_ref,V_rest;
  double ***trj;
  double  temp;

  char *inputfilename1,*inputfilename2,*parmfilename;
  char *outputfilename,*outputfilename2;

  FILE	*inputfile1,*inputfile2,*outputfile,*outputfile2,*parmfile;

  if (argc < 6) {
    printf("USAGE: %s [at] nwindow nt V_rest temp inputfilename2(parm) outputfilename(free energy) outputfilename2(const energy ) inputfilenames(trj) inputfilenames(crd ref)\n",argv[0]);
    exit(1);
  }
  flag=(*++argv)[0];
  if (flag != 'a' && flag != 't') {
    printf("flag error: must be   or  ");
    exit(1);
  }
  if (flag=='a') AMBERFLAG=ON;
  nwindows=atoi(*++argv);
  nttemp=atoi(*++argv);
  nt=(int *)gcemalloc(sizeof(int)*nwindows);
  for (i=0;i<nwindows;++i)
    nt[i]=nttemp;
  V_rest=atof(*++argv);
  temp=atof(*++argv);

  parmfilename=*++argv;
  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;
  numca=countatomtype("CA",2);

  outputfilename=*++argv;
  outputfilename2=*++argv;

  trj=(double ***)gcemalloc(sizeof(double **)*nwindows);
  for (i=0;i<nwindows;++i) {
    inputfilename1=*++argv;
    inputfile1=efopen(inputfilename1,"r");
    trj[i]=(double **)gcemalloc(sizeof(double *)*nt[i]);
    for (j=0;j<nt[i];++j)
      trj[i][j]=(double *)gcemalloc(sizeof(double)*numca*3);
    if (AMBERFLAG==ON) { getline(&line,&len,inputfile1); }
    io_scancatrj(inputfile1,numatom,nt[i],trj[i]);
    fclose(inputfile1);
  }

  crd_ref=(double **)gcemalloc(sizeof(double *)*nwindows);
  for (i=0;i<nwindows;++i)
    crd_ref[i]=(double *)gcemalloc(sizeof(double)*numca*3);
  for (i=0;i<nwindows;++i) {
    inputfilename2=*++argv;
    inputfile2=efopen(inputfilename2,"r");
    if (AMBERFLAG==ON) { getline(&line,&len,inputfile2);getline(&line,&len,inputfile2);  }
    io_scancaconf(inputfile2,numatom,crd_ref[i]);
    fclose(inputfile2);
  }

  ebF  = (double *)gcemalloc(sizeof(double)*nwindows);
    
  ebW  = (double ***)gcemalloc(sizeof(double **)*nwindows);
  for (i=0;i<nwindows;++i) {
    ebW[i]=(double **)gcemalloc(sizeof(double *)*nwindows);
    for (j=0;j<nwindows;++j)
      ebW[i][j]=(double *)gcemalloc(sizeof(double)*nt[j]);
  }
  
  wham_calc_force_carest(nwindows,nt,ebW,trj,numca,crd_ref,V_rest,temp);
  wham_ite(nwindows,nt,ebF,ebW);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<nwindows;++i)
    fprintf(outputfile,"%e\n",ebF[i]);  
  fclose(outputfile);

  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<nwindows;++i) {
    for (j=0;j<nwindows;++j) {
      for (k=0;k<nt[j];++k) {
	fprintf(outputfile2,"%e \n",ebW[i][j][k]);
      }
      fprintf(outputfile2,"\n");
    }
  }
  fclose(outputfile2);
 
  return 0;
}

