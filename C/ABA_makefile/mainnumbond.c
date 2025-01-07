#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "PT.h"

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i;
  double nb;

  int IHflag=EXCH;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *parmfilename,*progname;
  FILE *parmfile;

  while((c=getopt(argc,argv,"hH"))!=-1) {
    switch(c) {
    case 'H':
      IHflag=INCH;
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
  parmfilename = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);

  nb=numbonds(IHflag);

  printf("%d\n",nb);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-H] INCLUDE Hatom flag\n");
  printf("[-h] help \n");
  printf("%s [-H] [-h] parmfilename\n",progname);
}
