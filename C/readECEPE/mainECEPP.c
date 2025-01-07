#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ECEPE.h"

#include "EF.h"
//#include "IO.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;

  struct ECEPE_parms ECEPE_p;
  struct pnb nb_p;

  char *progname;
  char *preofilename,*bd8filename;

  char *line,*dummy;
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
  preofilename  = *argv;
  bd8filename = *++argv;

  read_ECEPE_parm(preofilename,bd8filename,&ECEPE_p,&nb_p);

  printf("NATOM=%4d ",ECEPE_p.NUMATM);
  printf("NDIHE=%4d ",ECEPE_p.NUMVAR);
  printf("NRESD=%4d ",ECEPE_p.NUMRES);
  printf("NINTE=%4d ",ECEPE_p.NUMINT);
  printf("NSS  =4%d\n",ECEPE_p.NUMS);

  for (i=0;i<ECEPE_p.NUMVAR;++i) {
    printf("A=%10.4lf ",ECEPE_p.dihed[i].A);
    printf("NB=%4d ",ECEPE_p.dihed[i].NB);
    printf("NS=%4d\n",ECEPE_p.dihed[i].NS);
  }

  for (i=0;i<ECEPE_p.NUMATM;++i) {
    for (j=0;j<3;++j)
      printf("%10.6lf ",ECEPE_p.atom[i].refcoord[j]);
    printf("charge=%10.6lf ",ECEPE_p.atom[i].charge);
    printf("nbtype=%d4 \n",ECEPE_p.atom[i].nbtype);
  }
  
  return 0;
}

void USAGE(char *progname) {
  printf("-h -- help\n");
  printf("USAGE: %s profilename bd8filename \n", progname);
}



