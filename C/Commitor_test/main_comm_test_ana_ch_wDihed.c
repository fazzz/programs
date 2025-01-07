#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <getopt.h>

#include "EF.h"

#define ON 1
#define OFF 0

int usage(void);

int main(int argc, char *argv[]) {
  int i,j,k,l,d,index,dummy;
  int numdihed=2,numspring;

  double *dihed_target,width;
  double dist,ddist,minddist;

  double *theta,**mfep;

  double pi;

  char *line;
  size_t len=0;

  extern char *optarg;
  extern int optind,opterr,optopt;

  int c;

  char *inputfilename;
  char *progname;
  FILE *inputfile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"numd",1,NULL,'d'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hd:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'd':
      numdihed=atoi(optarg);
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 2+numdihed) {
    USAGE(progname);
    exit(1);
  }
  numspring=atoi(*argv);
  theta=(double *)gcemalloc(sizeof(double)*numdihed);
  for (i=0;i<numdihed;++i) {
    theta[i]=atof(*++argv);
  }
  inputfilename   = *++argv;

  pi=acos(-1.0);

  inputfile=efopen(inputfilename,"r");
  mfep=(double **)gcemalloc(sizeof(double *)*numspring);
  for (i=0;i<numspring;++i) {
    mfep[i]=(double *)gcemalloc(sizeof(double)*numdihed);
  }
  for (i=0;i<numspring;++i) {
    for (j=0;j<numdihed;++j) {
      fscanf(inputfile,"%lf",&mfep[i][j]);
    }
  }
  fclose(inputfile);

  for (i=0;i<numspring;++i) {
    ddist=0.0;
    for (j=0;j<numdihed;++j) {
      dist=(theta[j]-mfep[i][j]);
      if (dist < 0.0) dist=-1.0*dist;
      if (dist > pi) dist=2.0*pi-dist;
      ddist+=dist*dist;
      ddist=sqrt(ddist);
    }
    if (i==0) {
      index=i;
      minddist=ddist;
    }
    else {
      if (minddist>ddist) {
	index=i;
	minddist=ddist;
      }
    }
  }

  printf("%d \n",index+1);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] \n");
  printf("%s inputfilename outputfilename \n",progname);
}

