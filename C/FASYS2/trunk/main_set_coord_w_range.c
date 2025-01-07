
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
  int num;

  char *line;
  size_t len=0;

  char *progname;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double pi,temp;
  double psi,phi;
  double psimin,psimax,psiint=10.0,phimin,phimax,phiint=10.0;
  double q[5][3];

  char *outputfilename;
  FILE *outputfile;

  progname=argv[0];

  while((c=getopt(argc,argv,"hi:j:"))!=-1) {
    switch(c) {
    case 'i':
      psiint=atof(optarg);
      break;
    case 'j':
      phiint=atof(optarg);
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

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  psimin=-atof(*argv);
  psimax=-atof(*++argv);
  phimin=-atof(*++argv);
  phimax=-atof(*++argv);
  outputfilename=*++argv;

  if (phimin>180.0 || phimin <-180.0) {
    printf("error:phi must be [-180:180]");
    exit(1);
  }

  if (psimin>180.0 || psimin <-180.0) {
    printf("error:psi must be [-180:180]");
    exit(1);
  }

  if (phimax>180.0 || phimax <-180.0) {
    printf("error:phi must be [-180:180]");
    exit(1);
  }

  if (psimax>180.0 || psimax <-180.0) {
    printf("error:psi must be [-180:180]");
    exit(1);
  }

  if (psimin > psimax) {
    temp=psimin;
    psimin=psimax;
    psimax=temp;
  }
  if (phimin > phimax) {
    temp=phimin;
    phimin=phimax;
    phimax=temp;
  }

  pi=acos(-1.0);

  psimin=psimin/180.0*pi;
  psimax=psimax/180.0*pi;
  psiint=psiint/180.0*pi;
  phimin=phimin/180.0*pi;
  phimax=phimax/180.0*pi;
  phiint=phiint/180.0*pi;

  num=0;
  outputfile=efopen(outputfilename,"w"); 
  for (phi=phimin;phi<phimax;phi+=phiint) {
    for (psi=psimin;psi<psimax;psi+=psiint) {
      ++num;
      set_coord_from_dihed(phi,psi,q);
      for (i=0;i<5;++i) {
	for (j=0;j<3;++j) 
	  fprintf(outputfile,"%4.1lf ",q[i][j]);
	fprintf(outputfile,"\n");
      }
      fprintf(outputfile,"\n");
    }
  }
  fprintf(outputfile,"%d\n",num);
  fclose(outputfile);
}

int USAGE(char *progname) {
  printf("USAGE:%s \n",progname);
  printf("-h \n");
  printf("-i psiint \n");
  printf("-j phiint \n");
  printf(" psimin psimax phimin phimax outputfilename\n");
}

 
