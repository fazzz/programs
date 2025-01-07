#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ECEPE.h"

#include "EF.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k;
  double f;

  double *dihed;

  struct ECEPE_parms ECEPE_p;
  struct pnb nb_p;

  char *progname;
  char *preofilename,*bd8filename;
  char *angfilename,*outputfilename;
  FILE *angfile,*outputfile;

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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  preofilename   = *argv;
  bd8filename    = *++argv;
  angfilename    = *++argv;
  outputfilename = *++argv;

  read_ECEPE_parm(preofilename,bd8filename,&ECEPE_p,&nb_p);

  dihed=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMVAR));

  angfile=efopen(angfilename,"r");
  fscanf(angfile,"%lf",&dihed[0]);
  for (i=0;i<9;++i) fscanf(angfile,"%lf",&f);
  for (i=1;i<ECEPE_p.NUMVAR;++i) fscanf(angfile,"%lf",&dihed[i]);
  fclose(angfile);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<ECEPE_p.NUMVAR;++i) fprintf(angfile,"%d %lf\n",i+1,dihed[i]);
  fclose(outputfile);

  return 0;
}

void USAGE(char *progname) {
  printf("-h -- help\n");
  printf("USAGE: %s profilename bd8filename angfilename outputfilename\n", progname);
}

