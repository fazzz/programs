
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PMF.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

int main(int argc, char *argv[]) {
  int i,j,flag=0;

  int numatom;
  double *data,*w;
  int numstep,numcol,xcol,ycol;
  double maxx,maxy,minx,miny;
  int framex,framey;
  double width,temp;
  double *pmf;
  
  char *inputfilename,*inputfilename2,*inputfilename3,*outputfilename,*c;
  FILE *inputfile1,*inputfile2,*inputfile3,*outputfile;
  
  if (argc < 5) {
    printf("USAGE: %s [at] temp width numstep numcol xcol ycol inputfilename2(parm) inputfilename(2d_data) inputfilename(w_func) outputfilename(pmf2d)   \n",argv[0]);
    exit(1);
  }
  flag=(*++argv)[0];
  if (flag != 'a' && flag != 't') {
    printf("flag error: must be a  or t ");
    exit(1);
  }
  temp=atof(*++argv);
  width=atof(*++argv);
  numstep=atoi(*++argv);
  numcol=atoi(*++argv);
  xcol=atoi(*++argv);
  ycol=atoi(*++argv);
  inputfilename2 = *++argv;
  inputfile2=efopen(inputfilename2,"r");
  readParmtop(inputfile2);
  fclose(inputfile2);
  inputfilename=*++argv;
  inputfile1=efopen(inputfilename,"r");
  data=(double *)gcemalloc(sizeof(double)*numstep*2);
  io_scnasdcoldata(inputfile1,numstep,0,numcol,xcol,ycol,data);
  fclose(inputfile1);
  w=(double *)gcemalloc(sizeof(double)*numstep);
  inputfilename3=*++argv;
  inputfile3=efopen(inputfilename3,"r");
  for (i=0;i<numstep;++i)
    fscanf(inputfile3,"%lf",&w[i]);
  fclose(inputfile3);
  numatom=AP.NATOM;
  outputfilename = *++argv;

  pmf=pmf_2dmap_rew(data,w,numstep,width,&maxx,&maxy,&minx,&miny,&framex,&framey,temp);
  
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<framex;++i) {
    for (j=0;j<framey;++j)
      fprintf(outputfile,"%e %e %e\n",width*i+minx,width*j+miny,pmf[i*framey+j]);
    fprintf(outputfile,"\n");
  }
  fclose(outputfile);
  
  return 0;
}

