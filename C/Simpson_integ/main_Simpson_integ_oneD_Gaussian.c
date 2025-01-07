
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "Simpson_integ.h"
#include "EF.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k;
  int N;

  double nyu=0.0,sigma=1.0;
  double S,minx,maxx;

  double pi;

  char *line;
  size_t len=0;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;
  
  pi=acos(-1.0);
  
  struct option long_opt[] = {
    {"nyu",1,NULL,'n'},
    {"sigma",1,NULL,'s'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"hn:s:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'n':
      nyu=atof(optarg);  break;
    case 's':
      sigma=atof(optarg);  break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);  exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  N = atoi(*argv);
  minx = atof(*++argv);
  maxx = atof(*++argv);

  S=Simpson_integ_oneD_Gaussian(N,minx,maxx,nyu,sigma,pi);

  printf("%10.8lf\n",S);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] [-nyu] average [-sigma] variance N minx miaxx \n",progname);
}
