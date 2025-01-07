#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "HIST.h"
#include "EF.h"

void USAGE(char *progname);

int scna1coldata(FILE *inputfile,int numf,int numi,int numcol,int col,double *data);

int main(int argc, char *argv[]) {
  int i,j;
  int numdata=1,numcol=1,col=1,normflag='o',numini=1;
  int frame;
  double width;
  double max,min;
  double *data,*hist;

  char *inputfilename,*outputfilename,*progname;
  FILE *inputfile,*outputfile;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  while((c=getopt(argc,argv,"hw:i:f:r:s:"))!=-1) {
    switch(c) {
    case 'n':
      normflag=optarg[0];
      break;
    case 'w':
      width=atof(optarg);
      break;
    case 'i':
      numini=atoi(optarg);
      break;
    case 'f':
      numdata=atoi(optarg);
      break;
    case 's':
      numcol=atoi(optarg);
      break;
    case 'r':
      col=atoi(optarg);
      break;
    case 'h':
      USAGE(progname);
      break;
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=argv[0];
  argc-=optind;
  argv+=optind;

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  inputfilename      =  *argv;
  outputfilename     =  *++argv;

  inputfile   =  efopen(inputfilename,"r");
  outputfile  =  efopen(outputfilename,"w");

  if (numdata-numini<=0) {
    printf("error: numdata is less than numini\n");
    exit(1);
  }

  data=(double *)gcemalloc(sizeof(double)*(numdata-numini+1));
  scna1coldata(inputfile,numdata,numini,numcol,col,data);
  hist=hist_mkhist(data,width,(numdata-numini+1),&max,&min,&frame,normflag);

  if (normflag=='f')
    for (i=0;i<=frame;++i) fprintf(outputfile,"%e %e \n",width*(i+0.5)+min,hist[i]);
  else if (normflag=='o')
    for (i=0;i<=frame;++i) fprintf(outputfile,"%e %e \n",width*(i+0.5)+min,hist[i]);
  else if (normflag=='a')
    for (i=0;i<=frame;++i) fprintf(outputfile,"%e %e \n",width*(i+0.5)+min,hist[i]/width);
  else {
    for (i=0;i<=frame;++i) {
      fprintf(outputfile,"%e %e \n",width*i+min,0.0);
      fprintf(outputfile,"%e %e \n",width*(i+1.0)+min,hist[i]/width);
    }
  }
  fclose(inputfile);
  fclose(outputfile);
}

int scna1coldata(FILE *inputfile,int numf,int numi,int numcol,int col,double *data){
  int i,j;
  double temp;
  char *line;
  size_t len=0;

  for (i=0;i<numi;++i)
    getline(&line,&len,inputfile);  
  for (i=0;i<numf-numi+1;++i) {
    for (j=0;j<numcol;++j) {
      if (j==col-1)
	fscanf(inputfile,"%lf",&data[i]);
      else 
	fscanf(inputfile,"%lf",&temp);
    }
  }

  return 0;
}

void USAGE(char *progname) {
  printf("[-n o/f/a/t ]  normflag");
  printf("[-w width]");
  printf("[-i numini]");
  printf("[-f numdata]");
  printf("[-s numcol]");
  printf("[-r col]");
  printf("%s inputfilename outputfilename\n",progname);
}
