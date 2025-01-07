
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "TOPO.h"
#include "LA.h"
#include "mymath.h"

#include "PTL.h"
#include "MB.h"
#include "EF.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numatom,na1,na2,numstep;
  
  double f,length,length_equ,V,E;
  double crd[2][3];

  double pi;

  char *inputfilename,*parmtopfilename,*outputfilename;
  FILE *inputfile,*parmfile,*outputfile;

  char *line;
  size_t len=0;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;
  
  pi=acos(-1.0);
  
  struct option long_opt[] = {
    {"na1",0,NULL,'h'},
    {"na2",0,NULL,'h'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"h",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);   exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < 8) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  parmtopfilename = *++argv;
  numstep = atoi(*++argv);
  na1 = atoi(*++argv);
  na2 = atoi(*++argv);
  length_equ = atof(*++argv);
  V = atof(*++argv);
  outputfilename = *++argv;

  parmfile=efopen(parmtopfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);

  outputfile=efopen(outputfilename,"w");

  for (i=0;i<numstep;++i) {
    for (j=0;j<numatom;++j) {
      if (j==na1-1) for (k=0;k<3;++k) fscanf(inputfile,"%lf",&crd[0][k]);
      else if (j==na2-1) for (k=0;k<3;++k) fscanf(inputfile,"%lf",&crd[1][k]);
      else for (k=0;k<3;++k) fscanf(inputfile,"%lf",&f);
    }

    length=0.0;
    for (j=0;j<3;++j) {
      length+=(crd[0][j]-crd[1][j])*(crd[0][j]-crd[1][j]);
    }
    length=sqrt(length);

    E=0.5*V*(length-length_equ)*(length-length_equ);
    fprintf(outputfile,"%d %10.4lf \n",i,E);

  }

  fclose(inputfile);
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parmtopfilename numstep na1 na2 length_equ V outputfilename \n",progname);
}
