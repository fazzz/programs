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
  int kflag=ON,jflag=OFF;

  int numatom;

  struct ECEPE_parms ECEPE_p;
  struct pnb nb_p;
  struct ECEPE_pote p;
  struct ECEPE_force f;

  char *name_atom,*name_res;

  char *progname;
  char *preofilename;
  FILE *parmtopfile;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"kjh"))!=-1) {
    switch(c) {
    case 'k':
      kflag=ON;
      break;
    case 'j':
      jflag=ON;
      break;
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
  numatom=ECEPE_p.NUMATM;

  name_atom=(char *)gcemalloc(sizeof(char)*numatom*4);
  name_res=(char *)gcemalloc(sizeof(char)*numatom*3);

  if (kflag==ON) {
    for (i=0;i<numatom;++i) {      
      for (j=0;j<4;++j)
	name_atom[(ECEPE_p.atom[i].katom-1)*4+j]=ECEPE_p.atom[i].name_atom[j];
      for (j=0;j<3;++j)
	name_res[(ECEPE_p.atom[i].katom-1)*3+j]=ECEPE_p.atom[i].name_res[j];
    }
  }


  for (i=0;i<numatom;++i) {
    if (kflag==ON) {
      for (j=0;j<4;++j) {
	printf("%c",name_atom[i*4+j]);
      }
      printf(" ");
      for (j=0;j<3;++j) {
	printf("%c",name_res[i*3+j]);
      }
      printf("\n");
    }
    else {
      printf("%s\n",ECEPE_p.atom[i].name_atom);
    }
  }

  return 0;
}

void USAGE(char *progname) {
  printf("-k -- k index\n");
  printf("-j -- j index\n");
  printf("-h -- help\n");
  printf("USAGE: %s profilename \n", progname);
}

