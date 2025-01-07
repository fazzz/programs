
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
  double a,b;
  double psi1,psi2,psiint=10.0,phi1,phi2;
  double q[5][3];

  char *outputfilename;
  FILE *outputfile;

  progname=argv[0];

  while((c=getopt(argc,argv,"hi:"))!=-1) {
    switch(c) {
    case 'i':
      psiint=atof(optarg);
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
  psi1=-atof(*argv);
  phi1=-atof(*++argv);
  psi2=-atof(*++argv);
  phi2=-atof(*++argv);
  outputfilename=*++argv;

  if (phi1 > 180.0 || phi1 < -180.0) {
    printf("error:phi must be [-180:180]");
    exit(1);
  }

  if (psi1 > 180.0 || psi1 < -180.0) {
    printf("error:psi must be [-180:180]");
    exit(1);
  }

  if (phi2 > 180.0 || phi2 < -180.0) {
    printf("error:phi must be [-180:180]");
    exit(1);
  }

  if (psi2 > 180.0 || psi2 < -180.0) {
    printf("error:psi must be [-180:180]");
    exit(1);
  }


  if (psi1 > psi2) {
    temp=psi1;
    psi1=psi2;
    psi2=temp;
  }

  pi=acos(-1.0);

  psi1=psi1/180.0*pi;
  psi2=psi2/180.0*pi;
  phi1=phi1/180.0*pi;
  phi2=phi2/180.0*pi;
  psiint=psiint/180.0*pi;

  a=(phi2-phi1)/(psi2-psi1);
  b=phi1-a*psi1;

  num=0;
  outputfile=efopen(outputfilename,"w"); 
  for (psi=psi1;psi<psi2;psi+=psiint) {
    phi=a*psi+b;
    ++num;
    set_coord_from_dihed(phi,psi,q);
    for (i=0;i<5;++i) {
      for (j=0;j<3;++j) 
	fprintf(outputfile,"%4.1lf ",q[i][j]);
      fprintf(outputfile,"\n");
    }
    fprintf(outputfile,"\n");
  }
  fprintf(outputfile,"%d\n",num);
  fclose(outputfile);
}

int USAGE(char *progname) {
  printf("USAGE:%s \n",progname);
  printf("-h \n");
  printf("-i psiint \n");
  printf(" psi1 phi1 psi1 psi2 outputfilename\n");
}

 
