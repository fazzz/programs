
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
  
  int nT,*T;
  double Tmin,Tmax,xmin,xmax,*x,dx;

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
  nT        = atoi(*++argv);

  T=(int *)gcemalloc(sizeof(int)*nT);
  x=(double *)gcemalloc(sizeof(double)*nT);

  xmin=log(Tmin);
  xmax=log(Tmax);
  dx=(xmax-xmin)/(nT-1);

  for (i=0;i<nT;++i) x[i]=xmin+dx*i;
  for (i=0;i<nT;++i) T[i]=(int)exp(x[i]);
  
  for (i=0;i<nT;++i) printf("%4d ",T[i]);
  printf("\n");

  return 0.0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] Tmin Tmax nT \n",progname);
}
