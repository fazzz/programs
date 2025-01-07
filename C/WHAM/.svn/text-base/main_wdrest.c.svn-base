
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

  int nwindows,numatom,*nt;
  int nttemp;
  double *ebF,***ebW;
  int **dihedpairs;
  int numdrest;
  double **theta_ref,V_rest;
  double ***trj;
  double  temp;

  char	 *inputfilename1,*inputfilename2,*parmfilename;
  char	 *outputfilename,*outputfilename2;

  FILE	*inputfile1,*inputfile2,*outputfile,*outputfile2,*parmfile;

  if (argc < 6) {
    printf("USAGE: %s nwindow nt numdrest V_rest temp inputfilename2(parm) outputfilename(fk) inputfilenames(trj) inputfilename2(parm_rest) \n",argv[0]);
    exit(1);
  }
  nwindows=atoi(*++argv);
  nttemp=atoi(*++argv);
  nt=(int *)gcemalloc(sizeof(int)*nwindows);
  for (i=0;i<nwindows;++i)
    nt[i]=nttemp;
  numdrest=atoi(*++argv);
  V_rest=atof(*++argv);
  temp=atof(*++argv);
  parmfilename=*++argv;
  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;
  outputfilename=*++argv;
  trj=(double ***)gcemalloc(sizeof(double **)*nwindows);
  for (i=0;i<nwindows;++i) {
    trj[i]=(double **)gcemalloc(sizeof(double *)*nt[i]);
    for (j=0;j<nt[i];++j)
      trj[i][j]=(double *)gcemalloc(sizeof(double)*numatom*3);
  }
  for (i=0;i<nwindows;++i) {
    inputfilename1=*++argv;
    inputfile1=efopen(inputfilename1,"r");
    io_scantrj(inputfile1,numatom,nt[i],trj[i]);
    fclose(inputfile1);
  }
  dihedpairs=(int **)gcemalloc(sizeof(int *)*nwindows);
  for (i=0;i<nwindows;++i)
    dihedpairs[i]=(int *)gcemalloc(sizeof(int)*numdrest*4);
  theta_ref=(double **)gcemalloc(sizeof(double *)*nwindows);
  for (i=0;i<nwindows;++i)
    theta_ref[i]=(double *)gcemalloc(sizeof(double)*numdrest);
  for (i=0;i<nwindows;++i) {
    //  fscanf(inputfile3,"%d",&numdrest); 
    inputfilename2=*++argv;
    inputfile2=efopen(inputfilename2,"r");
    for (j=0;j<numdrest;++j)
      for (k=0;k<4;++k)
	fscanf(inputfile2,"%d",&dihedpairs[i][j*4+k]);
    for (j=0;j<numdrest;++j)
      fscanf(inputfile2,"%lf",&theta_ref[i][j]);
    fclose(inputfile2);
  }

  ebF  = (double *)gcemalloc(sizeof(double)*nwindows);
    
  ebW  = (double ***)gcemalloc(sizeof(double **)*nwindows);
  for (i=0;i<nwindows;++i) {
    ebW[i]=(double **)gcemalloc(sizeof(double *)*nwindows);
    for (j=0;j<nwindows;++j) {
      ebW[i][j]=(double *)gcemalloc(sizeof(double)*nt[j]);
    }
  }
  
  wham_calc_force_drest(nwindows,nt,ebW,trj,numatom,dihedpairs,theta_ref,numdrest,V_rest,temp);
  wham_ite(nwindows,nt,ebF,ebW);
  //  wham_pmf(num_k,num_x,num_y,num_window,nt,ebF,ebW,hist,pmf);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<nwindows;++i)
    fprintf(outputfile,"%e\n",ebF[i]);  
  fclose(outputfile);

  outputfile2=efopen("log_ef.txt","w");
  for (i=0;i<nwindows;++i) {
    for (j=0;j<nwindows;++j) {
      for (k=0;k<nt[j];++k) {
	fprintf(outputfile2,"%e ",ebW[i][j][k]);
      }
      fprintf(outputfile2,"\n");
    }
  }
  fclose(outputfile2);

 
  return 0;
}

