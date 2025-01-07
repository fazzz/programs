#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "HIST.h"
#include "EF.h"

int scnadcoldata(FILE *inputfile,int numf,int numi,int numcol,int xcol,int ycol, double *data);

int main(int argc, char *argv[]) {
  int i,j;
  int numdata,numcol,xcol,ycol,normflag,logflag;
  int framex,framey;
  double width;
  double maxx,maxy,minx,miny;
  double *data,*hist;

  char *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;

  if (argc < 7) {
    printf("USAGE: ./%s normflag( o or f  ) inputfilename(data) outputfilename(hist) width numdata numcol xcol ycol\n",argv[0]);
    exit(1);
  }
  normflag=(*++argv)[0];
  if (normflag != 'o' && normflag != 'f') {
    printf("flag error: must be o(n)  or (of)f ");
    exit(1);
  }
  else if (normflag == 'o')
    normflag=ON;
  else
    normflag=OFF;
  inputfilename   = *++argv;
  outputfilename  = *++argv;
  width = atof(*++argv);
  numdata =atoi(*++argv);
  numcol =atoi(*++argv);
  if (argc < 7)
    xcol=1;
  else
    xcol = atoi(*++argv);
  if (argc < 8)
    ycol = atoi(*++argv);
  else
    ycol = atoi(*++argv);

  inputfile   =  efopen(inputfilename,"r");
  outputfile  = efopen(outputfilename,"w");

  data=(double *)gcemalloc(sizeof(double)*numdata*2);

  scnadcoldata(inputfile,numdata,1,numcol,xcol,ycol,data);
  hist=hist_mk2dhist(data,width,width,numdata,&maxx,&maxy,&minx,&miny,&framex,&framey,normflag);

  for (i=0;i<=framex;++i) {
    for (j=0;j<framey;++j)
      fprintf(outputfile,"%e %e %e\n",width*i+minx,width*j+miny,hist[i*framey+j]);
    fprintf(outputfile,"\n");
  }

  fclose(inputfile);
  fclose(outputfile);
}

int scnadcoldata(FILE *inputfile,int numf,int numi,int numcol,int xcol,int ycol, double *data){
  int i,j;
  double temp;
  char *line;
  size_t len=0;

  for (i=0;i<numi;++i)
    getline(&line,&len,inputfile);  
  for (i=0;i<numf-numi+1;++i) {
    for (j=0;j<numcol;++j) {
      if (j==xcol-1)
	fscanf(inputfile,"%lf",&data[i*2]);
      else if (j==ycol-1)
	fscanf(inputfile,"%lf",&data[i*2+1]);
      else 
	fscanf(inputfile,"%lf",&temp);
    }
  }

  return 0;
}

