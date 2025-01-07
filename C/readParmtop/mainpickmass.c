#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "PT.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j;
  int flag=OFF;
  int numatom;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *parmfilename,*progname;
  FILE *parmfile;

  while((c=getopt(argc,argv,"hn"))!=-1) {
    switch(c) {
    case 'n':
      flag=ON;
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

  numatom=AP.NATOM;

  j=0;
  for (i=0;i<numatom;++i) {
    if (flag==ON) {
      if ((i)==AP.IPRES[j]-1 ) {
	++j;
	printf("\n");
      }
    }
    else if ((i+1)%10==0) printf("\n");
    printf("%7.3lf ",AP.AMASS[i]);
  }
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-n] flag on \n");
  printf("[-h] help \n");
  printf("%s [-h] parmfilename\n",progname);
}
