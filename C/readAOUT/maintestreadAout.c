#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "readAOUT.h"
#include "EF.h"

int usage(char *progname);

int main(int argc, char *argv[]) {
  int i;
  int numstep;
  double *data;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*outputfilename;
  char *progname;
  FILE *inputfile,*outputfile;

  while((c=getopt(argc,argv,"h"))!=-1) {
    switch(c) {
    case 'h':
      usage(progname);
      exit(1);
    default:
      usage(progname);
      exit(1);
    }
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 2) {
    usage(progname);
    exit(1);
  }
  inputfilename = *argv;
  outputfilename = *++argv;

  inputfile=efopen(inputfilename,"r");
  data=readOUT(inputfile,"toal_vertial_energy      =",&numstep);
  fclose(inputfile);
  
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    fprintf(outputfile,"%8.5e\n",data[i]);
  }
  fclose(outputfile);
}

int usage(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s inputfilename  outputfilename \n",progname);
}
