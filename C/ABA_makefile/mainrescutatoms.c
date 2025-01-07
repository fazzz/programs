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
  int i,j,d;
  int numatom,ires;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crd;

  char *inputfilename;
  char *parmfilename,*outputfilename,*progname;
  FILE *parmfile,*outputfile,*inputfile;

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

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  parmfilename = *argv;
  inputfilename = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);

  numatom=AP.NATOM;

  inputfile=efopen(inputfilename,"r");
  crd=(double *)gcemalloc(sizeof(double )*numatom*3);
  fscanf(inputfile,"%d",&d);
  ires=0;
  for (i=0;i<numatom;++i) {
    if (AP.IPRES[ires+1]<i) {
      ires+=1;
      printf("\n"); 
    }
    for (j=0;j<3;++j) {
      fscanf(inputfile,"%lf",&crd[i*3+j]); 
      printf("%lf ",crd[i*3+j]); 
    }
    for (j=0;j<4;++j) {
      printf("%c ",AP.IGRAPH[i][j]); 
    }
    printf(" "); 
    for (j=0;j<3;++j) {
      printf("%c",AP.LABERES[ires][j]); 
    }
    printf("\n"); 
  }

  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] parmfilename \n",progname);
}
