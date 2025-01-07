
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k;
  
  int M;
  double Tmin,Tmax,f;

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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  Tmin      = atof(*argv);
  Tmax      = atof(*++argv);
  f         = atof(*++argv);

  M = (int)(sqrt(f)*log(Tmax/Tmin));

  printf("%4d \n",M);

  return 0.0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] Tmin Tmax f (degree of freedom) \n",progname);
}
