#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ECEPE.h"
#include "PT.h"
#include "FF.h"

#include "EF.h"

#define ON  1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k;
  int dflag=OFF,aflag=OFF,nflag=OFF;

  int numdih,numatom;

  struct ECEPE_parms ECEPE_p;
  struct pnb nb_p;
  struct ECEPE_pote p;
  struct ECEPE_force f;

  char *progname;
  char *preofilename;
  FILE *parmtopfile;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"hdan"))!=-1) {
    switch(c) {
    case 'd':
      dflag=ON;
      break;
    case 'a':
      aflag=ON;
      exit(1);
    case 'n':
      nflag=ON;
      exit(1);
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

  if (argc < 1) {
    USAGE(progname);
    exit(1);
  }
  preofilename  = *argv;

  read_ECEPE_parm_wobd8(preofilename,&ECEPE_p);

  if ( aflag == OFF && dflag == OFF ) {
    USAGE(progname);
    exit(1);
  }
  else if ( aflag == ON  && dflag == ON )
    printf("%d %d",ECEPE_p.NUMATM,ECEPE_p.NUMVAR);
  else if ( aflag == ON )
    printf("%d",ECEPE_p.NUMATM);
  else 
    printf("%d",ECEPE_p.NUMVAR);
  if (nflag==OFF)
    printf("\n");

  return 0;
}

void USAGE(char *progname) {
  printf("-d -- numdih\n");
  printf("-a -- numatom\n");
  printf("-n -- wo new line\n");
  printf("-h -- help\n");
  printf("USAGE: %s profilename \n", progname);
}

