
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "SEQUENCE.h"
#include "EF.h"

#define ON 0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k;
  
  int num=10,nseq;
  double min=0.0,max=10.0,common=1.2;
  double *seq;

  int maxflag=ON;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"min",1,NULL,'i'},
    {"max",1,NULL,'a'},
    {"num",1,NULL,'n'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hi:a:n:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'i':
      min=atof(optarg);
      break;
    case 'a':
      max=atof(optarg);
      break;
    case 'n':
      num=atoi(optarg);
      maxflag=OFF;
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

  if (argc < 1) {
    USAGE(progname);
    exit(1);
  }
  common    = atof(*argv);

  if (common<=0.0) {
    printf("error\n");
  }

  if (maxflag==ON){
    if (max<=min)    printf("error\n");
  }
  else {
    if (num<=0)  printf("error\n");
  }

  if (maxflag==ON) seq=geome_seq_wmin_max(min,max,common);
  else seq=geome_seq_wmin_num(min,num,common);

  if (maxflag==ON)
    nseq=sizeof(seq)/sizeof(double);
  else 
    nseq=num;
  
  for (i=0;i<nseq;++i) printf("%4.2lf ",seq[i]);  printf("\n");

  return 0.0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] Tmin Tmax nT \n",progname);
}
