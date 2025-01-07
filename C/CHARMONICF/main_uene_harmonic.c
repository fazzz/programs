
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EF.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k;
  int numstep;

  double Kx, Ky, x0, y0;
  double *x,*y,p;

  double dx,dy;
  double pi;

  char  *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;
  
  char *line;
  size_t len=0;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;
  
  struct option long_opt[] = {
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"h",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);   exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < 7) {
    USAGE(progname);
    exit(1);
  }
  numstep = atoi(*argv);
  Kx = atof(*++argv);
  Ky = atof(*++argv);
  x0 = atof(*++argv);
  y0 = atof(*++argv);
  inputfilename = *++argv;
  outputfilename = *++argv;

  pi=acos(-1.0);

  x=(double *)gcemalloc(sizeof(double)*numstep);
  y=(double *)gcemalloc(sizeof(double)*numstep);

  inputfile=efopen(inputfilename,"r");
  for (i=0;i<numstep;++i) {
    fscanf(inputfile,"%lf %lf",&x[i],&y[i]);
  }
  fclose(inputfile);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {

    if ((dx=x[i]-x0)>pi)
      dx-=2.0*pi;
    else if ((dx=x[i]-x0)<-1.0*pi)
      dx+=2.0*pi;

    if ((dy=y[i]-y0)>pi)
      dy-=2.0*pi;
    else if ((dy=y[i]-y0)<-1.0*pi)
      dy+=2.0*pi;

    p=0.5*k*dx*dx+0.5*k*dy*dy;
    fprintf(outputfile,"%8.3lf\n",p);
  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] numstep Kx Ky x0 y0 inputfilename outputfilename\n",progname);
}


