
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "FASYS.h"
#include "EF.h"
#include "IO.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int flag;

  char *line;
  size_t len=0;

  char *progname;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double pi;
  double psi,phi,q[5][3];

  char *outputfilename;
  FILE *outputfile;

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
  phi=atof(*argv);
  psi=atof(*++argv);
  outputfilename=*++argv;

  if (phi>180.0 || phi <-180.0) {
    printf("error:phi must be [-180:180]");
  }

  if (psi>180.0 || psi <-180.0) {
    printf("error:psi must be [-180:180]");
  }

  pi=acos(-1.0);
  phi=phi/180.0*pi;
  psi=psi/180.0*pi;

  outputfile=efopen(outputfilename,"w"); 
  set_coord_from_dihed(phi,psi,q);
  for (i=0;i<5;++i) {
    for (j=0;j<3;++j) 
      fprintf(outputfile,"%10.4lf",q[i][j]);
    fprintf(outputfile,"\n");
  }
  fclose(outputfile);
}

int USAGE(char *progname) {
  printf("USAGE:%s \n",progname);
  printf(" psi phi\n");
}

 
