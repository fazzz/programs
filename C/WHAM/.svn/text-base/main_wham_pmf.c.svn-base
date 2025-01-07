
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "WHAM.h"
#include "HIST.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

int main(int argc, char *argv[]) {
  int i,j,k;

  double width,maxx,maxy,minx,miny,minx2,miny2;
  int framex,framey,framex2,framey2;
  int num_k;
  int numcol,xcol,ycol;
  
  int nwindows,numatom,*nt;
  int nttemp;
  double *ebF,***ebW;
  double **hist_ts;
  double **data;
  double **pmf;

  char *inputfilename1,*inputfilename2,*inputfilename3,*inputfilename4;
  char *outputfilename;
  FILE	*inputfile1,*inputfile2,*inputfile3,*inputfile4,*outputfile;

  if (argc < 6) {
    printf("USAGE: %s nwindow nt width num_k numcol xcol ycol inputfilename2(parm) inputfilename3(ebW) inputfilename4(ebF) outputfilename(pmf2d) inputfilenames(2ddata) \n",argv[0]);
    exit(1);
  }
  nwindows=atoi(*++argv);
  nttemp=atoi(*++argv);
  nt=(int *)gcemalloc(sizeof(int)*nwindows);
  for (i=0;i<nwindows;++i) nt[i]=nttemp;
  width=atof(*++argv);
  num_k=atoi(*++argv);
  numcol=atoi(*++argv);
  xcol=atoi(*++argv);
  ycol=atoi(*++argv);

  inputfilename2=*++argv;
  inputfile2=efopen(inputfilename2,"r");
  readParmtop(inputfile2);
  fclose(inputfile2);
  numatom=AP.NATOM;

  ebF  = (double *)gcemalloc(sizeof(double)*nwindows);
  ebW  = (double ***)gcemalloc(sizeof(double **)*nwindows);
  for (i=0;i<nwindows;++i) {
    ebW[i]=(double **)gcemalloc(sizeof(double *)*nwindows);
    for (j=0;j<nwindows;++j)
      ebW[i][j]=(double *)gcemalloc(sizeof(double)*nt[j]);
  }

  inputfilename3=*++argv;
  inputfile3=efopen(inputfilename3,"r");
  for (i=0;i<nwindows;++i)
    for (j=0;j<nwindows;++j)
      for (k=0;k<nt[j];++k)
	fscanf(inputfile3,"%lf",&ebW[i][j][k]);
  fclose(inputfile3);

  inputfilename4=*++argv;
  inputfile4=efopen(inputfilename4,"r");
  for (i=0;i<nwindows;++i)
    fscanf(inputfile4,"%lf",&ebF[i]);
  fclose(inputfile4);

  outputfilename=*++argv;
  data=(double **)gcemalloc(sizeof(double *)*nwindows);
  for (i=0;i<nwindows;++i)
    data[i]=(double *)gcemalloc(sizeof(double)*nt[i]*2);
  for (i=0;i<nwindows;++i) {
    inputfilename1=*++argv;
    inputfile1=efopen(inputfilename1,"r");
    io_scnadcoldata(inputfile1,nt[i],2,numcol,xcol,ycol,data[i]);
    fclose(inputfile1);
  }  

  hist_ts=(double **)gcemalloc(sizeof(double *)*nwindows);
  framex=0;
  framey=0;
  for (i=0;i<nwindows;++i) {
    hist_ts[i]=hist_mk2dhist_ts(data[i],width,width,nt[i],&maxx,&maxy,&minx2,&miny2,&framex2,&framey2,1);
    if (framex < framex2) framex=framex2;
    if (framey < framey2) framey=framey2;
    if (minx > minx2 || i==0) minx=minx2;
    if (miny > miny2 || i==0) miny=miny2;
  }

  pmf=(double **)gcemalloc(sizeof(double *)*framex);
  for (i=0;i<framex;++i) pmf[i]=(double *)gcemalloc(sizeof(double)*framey);
  wham_pmf2d(num_k-1,framex,framey,nwindows,nt,ebF,ebW,hist_ts,pmf);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<framex;++i) {
    for (j=0;j<framey;++j)
      fprintf(outputfile,"%e %e %e\n",width*i+minx,width*j+miny,pmf[i][j]);
    fprintf(outputfile,"\n");
  }
  fclose(outputfile);

  return 0;
}
