
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <getopt.h>

#include "EF.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i;

  int interval=10;

  int c;
  char *line;
  size_t len=0;

  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"int",1,NULL,'i'},
    {"h",1,NULL,'h'},
    {0,0,0,0}
  };

  progname=*argv;

  while((c=getopt_long(argc,argv,"hi:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'i':
      interval=atoi(optarg);
      break;
    default:
      USAGE(progname);
      exit(1);
    }
  }

  argc -= optind;
  argv += optind;

  if (argc < 2 ) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  outputfilename = *++argv;
  
  inputfile=efopen(inputfilename,"r");
  outputfile=efopen(outputfilename,"w");
  i=0;
  while (getline(&line, &len, inputfile) != -1) {
    if (i%interval==0) {
      fprintf(outputfile,"%s", line);
    }
    ++i;
  }

  fclose(inputfile);
  fclose(outputfile);

  return 0;
}

void USAGE(char *progname){
  printf("[-h] \n");
  printf("%s [-h] inputfilename outputfilename\n",progname);
}
