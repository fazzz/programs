#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "FFL.h"
#include "EF.h"

#define ON 1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i;
  int numstep=1000000;
  double omega=0.1,d=0.001;
  double pi;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *outputfilename;
  FILE *outputfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"d",1,NULL,'d'},
    {"nums",1,NULL,'s'},
    {"omega",1,NULL,'o'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hd:s:o:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'd':
      d=atof(optarg);
      break;
    case 'o':
      omega=atof(optarg);
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

  pi=acos(-1.0);

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 1) {
    USAGE(progname);
    exit(1);
  }
  outputfilename    = *argv;

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    v=0.0;
    for (j=0;j<num;++j) {
      v+=sin(2.0*pi/omega[j]*d*i)
    }
    fprintf(outputfile,"%e\n",v);
  }
  fclose(outputfile);

  return 0;
}

void USAGE(char *progname){
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] outputfilename\n",progname);
}

