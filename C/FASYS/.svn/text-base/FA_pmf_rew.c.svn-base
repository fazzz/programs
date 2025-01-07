
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PMF.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

double kbT=1.98723e-3*4.18407*100.0;

int main(int argc, char *argv[]) {
  int i,j,d;

  double beta;
  double *data,*w;
  int numstep,numcol,xcol,ycol;
  int numcol2,col2;
  double maxx,maxy,minx,miny;
  int framex,framey;
  double width,temp;
  double *pmf;
  
  char *datafilename,*enefilename,*outputfilename,*c;
  FILE *datafile,*enefile,*outputfile;
  
  if (argc < 5) {
    printf("USAGE: %s beta width numstep numcol xcol ycol numcol2 col2 inputfilename(2d_data) inputfilename(enedat) outputfilename(pmf2d)\n",argv[0]);
    exit(1);
  }
  beta=atof(*++argv);
  width=atof(*++argv);
  numstep=atof(*++argv);
  numcol=atoi(*++argv);
  xcol=atoi(*++argv);
  ycol=atoi(*++argv);
  numcol2=atoi(*++argv);
  col2=atoi(*++argv);
  datafilename = *++argv;
  enefilename  = *++argv;
  outputfilename = *++argv;

  datafile=efopen(datafilename,"r");
  data=(double *)gcemalloc(sizeof(double)*numstep*2);
  io_scnasdcoldata(datafile,numstep,0,numcol,xcol,ycol,data);
  fclose(datafile);

  w=(double *)gcemalloc(sizeof(double)*numstep);
  enefile=efopen(enefilename,"r");
  io_scancoldata(datafile,numstep,0,numcol2,col2,w);
  fclose(enefile);

  temp=1.0/kbT/beta;
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

