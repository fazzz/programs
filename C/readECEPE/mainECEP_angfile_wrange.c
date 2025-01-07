#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ECEPE.h"
#include "PROTOPO.h"
#include "PT.h"
#include "FF.h"
#include "TOPO.h"

#include "EF.h"

#define ON  1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k;

  double f;

  char *progname;
  char *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;

  char **line,*dummy;
  char **line_dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"h"))!=-1) {
    switch(c) {
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  argc-=optind;
  argv+=optind;

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  outputfilename = *++argv;

  inputfile=efopen(inputfilename,"r");
  outputfile=efopen(outputfilename,"w");
  while((c=fscanf(inputfile,"%lf",&f))!=-1) {
    if (f>180.0) f-=360.0;
    else if (f<-180.0) f+=360.0;
    fprintf(outputfile,"%lf\n",f);
  }
  fclose(inputfile);
  fclose(outputfile);

  return 0;
}

void USAGE(char *progname) {
  printf("-h -- help\n");
  printf("USAGE: %s inputfilename outputfilename \n", progname);
}

