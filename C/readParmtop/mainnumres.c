#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "PT.h"

void usage(char *progname);

#define ON 1
#define OFF 0

int main(int argc, char *argv[]) {
  int i;
  int nres;
  int nflag=ON;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *parmfilename,*progname;
  FILE *parmfile;

  while((c=getopt(argc,argv,"hn"))!=-1) {
    switch(c) {
   case 'n':
     nflag=OFF;
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
  parmfilename = *argv;

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);

  nres=AP.NRES;

  printf("%d",nres);
  if (nflag==ON)
    printf("\n",nres);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-n] new line \n");
  printf("[-h] help \n");
  printf("%s [-h] parmfilename\n",progname);
}
