#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "PTL.h"

void usage(char *progname);

#define ON 1
#define OFF 0

int main(int argc, char *argv[]) {
  int i,j;
  int natom,ires;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *parmfilename,*progname;
  FILE *parmfile;

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

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 1) {
    USAGE(progname);
    exit(1);
  }
  parmfilename = *argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile);

  natom=AP.NATOM;
  printf("%d\n",natom);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] parmfilename \n",progname);
}
