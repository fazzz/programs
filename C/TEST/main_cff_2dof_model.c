
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EF.h"

#define ON 1
#define OFF 0

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j;
  int numstep=100000;
  int addflag=OFF;
  
  double x,y;
  double V;

  double beta=1,kb=1.98723e-3;

  double D=5.0,a=1.0,K=1.0,lamda=2.878;

  char *outputfilename,*xtrjfilename,*ytrjfilename;
  FILE *outputfile,*xtrjfile,*ytrjfile;

  char *progname;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  int opt_idx=1;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  struct option long_opt[] = {
    {"temp",1,NULL,'t'},
    {"a",0,NULL,'a'},
    {"dx",1,NULL,'d'},
    {"int",1,NULL,'i'},
    {"nums",1,NULL,'s'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"has:d:i:t:r:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 't':
      beta=atof(optarg);
      break;
    case 's':
      numstep=atoi(optarg);
      break;
    case 'a':
      addflag=ON;
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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  xtrjfilename       = *argv;
  ytrjfilename       = *++argv;
  outputfilename     = *++argv;

  beta=1.0/(kb*beta);

  xtrjfile=efopen(xtrjfilename,"r");
  ytrjfile=efopen(ytrjfilename,"r");
  if (addflag==OFF)  outputfile=efopen(outputfilename,"w");
  if (addflag==ON)   outputfile=efopen(outputfilename,"a");
  for (i=0;i<numstep;++i) {
    fscanf(xtrjfile,"%lf ",&x);
    fscanf(ytrjfile,"%lf ",&y);

    V=D*(x*x-a*a)*(x*x-a*a)+0.5*K*y*y+lamda*x*y;
    // V=0.5*D*(x-a)*(x-a)+0.5*K*y*y;

    fprintf(outputfile,"%10.8lf \n",beta*V);
  }
  fclose(outputfile);
  fclose(xtrjfile);
  fclose(ytrjfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] -beta outputfilename trjfilename outrstfilename\n",progname);
}


