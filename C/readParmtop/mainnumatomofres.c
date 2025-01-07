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
  readParmtop(parmfile);
  fclose(parmfile);

  nres=AP.NRES;

  for (i=0;i<nres-1;++i) {
    if ((i+1)%10==0) printf("\n"); 
    printf("%d(%s) ",AP.IPRES[i+1]-AP.IPRES[i],AP.LABERES[i]);
  }
  if ((i+1)%10==0) printf("\n"); 
  printf("%d(%s) ",AP.NATOM-AP.IPRES[i],AP.LABERES[i]);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] parmfilename\n",progname);
}
