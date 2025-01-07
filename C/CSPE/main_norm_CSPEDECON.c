
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>

#include "fftw3.h"

#include "IO.h"
#include "PT.h"
#include "EF.h"
#include "SPE.h"
#include "const.h"

#include "netcdf_mine.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,n;
  double f;

  int numdihed=0;
  int numstep=0;

  double cv=2.999792e-2;
  double deltat=0.001,pi,kbT,temp=300;
  double *spe,*sumspe;

  char *inputfilename,*outputfilename,*progname;
  FILE *inputfile,*outputfile;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *line;
  size_t len=0;

  while((c=getopt(argc,argv,"hcvn:d:s:i:f:t:"))!=-1) {
    switch(c) {
    case 'n':
      numdihed=atoi(optarg);
      break;
    case 'd':
      deltat=atof(optarg);
      break;
    case 's':
      numstep=atoi(optarg);
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  if (numdihed==0) {
    printf("numdihed must not be 0\n");
    exit(1);
  }
  else if (numstep==0) {
    printf("num must not be 0\n");
    exit(1);
  } 

  progname=*argv;

  argc-=optind;
  argv+=optind;
  
  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  outputfilename = *++argv;

  pi=acos(-1.0);
  
  spe=(double *)gcemalloc(sizeof(double)*numstep*numdihed);
  
  inputfile=efopen(inputfilename,"r");
  for (i=0;i<numstep;++i) {
    fscanf(inputfile,"%lf",&f);
    for (j=0;j<numdihed;++j)
      fscanf(inputfile,"%lf",&spe[i*numdihed+j]);
  }
  fclose(inputfile);

  sumspe=(double *)gcemalloc(sizeof(double)*numdihed);
  for (i=0;i<numdihed;++i) sumspe[i]=0.0;
  for (i=0;(double)i/numstep/deltat/cv<1000.0;++i)
    for (j=0;j<numdihed;++j) sumspe[j]+=spe[i*numdihed+j];
  for (i=0;(double)i/numstep/deltat/cv<1000.0;++i)
    for (j=0;j<numdihed;++j) spe[i*numdihed+j]/=sumspe[j];

  outputfile=efopen(outputfilename,"w");
  pi=acos(-1.0);
  for (i=0;i<numstep;++i) {
    fprintf(outputfile,"%e ",(double)i/numstep/deltat/cv);
    for (j=0;j<numdihed;++j)
	fprintf(outputfile,"%e ",spe[i*numdihed+j]);
    fprintf(outputfile,"\n ");
  }
  fclose(outputfile);

  return 0;
}

void USAGE( char *progname) {
  printf("[-h] help\n");
  printf("[-n numdihed ]\n");
  printf("[-s numstep ]\n");
  printf("USAGE: %s inputfilename outputfilename\n",progname);
}
