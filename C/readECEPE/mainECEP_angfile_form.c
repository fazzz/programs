#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ECEPE.h"
#include "PROTOPO.h"
#include "PT.h"
#include "FF.h"
#include "TOPO.h"

#include "EF.h"

#define ON  1
#define OFF 0

#define MAXATOM 1000

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,m,n;

  double f;

  int PROflag=OFF;

  int *numreslist;

  char *progname;
  char *preofilename,*bd8filename,*inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;

  char **line,*dummy;
  char **line_dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  int ntotaldih,ndihinres,numres;

  struct ECEPE_parms ECEPE_p;
  struct pnb nb_p;


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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  preofilename   = *argv;
  bd8filename    = *++argv;
  inputfilename  = *++argv;
  outputfilename = *++argv;

  read_ECEPE_parm_wtransindex(preofilename,bd8filename,&ECEPE_p,&nb_p);

  ntotaldih=0;
  ndihinres=1;
  numres=1;
  numreslist=(int *)gcemalloc(sizeof(int)*ECEPE_p.NUMRES);
  for (i=0;i<ECEPE_p.NUMVAR;++i) {
    if (ECEPE_p.dihed[i].indexv1>numres) {
      numres=ECEPE_p.dihed[i].indexv1;
      numreslist[numres-2]=ndihinres;
      ndihinres=1;
    }
    if (ECEPE_p.dihed[i].indexv2>ndihinres) {
      ndihinres=ECEPE_p.dihed[i].indexv2;
    }
  }

  inputfile=efopen(inputfilename,"r");
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<ECEPE_p.NUMVAR;++i) {
    for (j=0;j<numreslist[i];++j) {
      fscanf(inputfile,"%lf",&f);
      fprintf(outputfile,"%lf\n",f);
    }
    for (j=0;j<10-numreslist[i];++j) {
      fscanf(inputfile,"%lf",&f);
    }
  }
  fclose(inputfile);
  fclose(outputfile);

  return 0;
}

void USAGE(char *progname) {
  printf("-h -- help\n");
  printf("USAGE: %s profilename bd8filename inputfilename outputfilename \n", progname);
}

