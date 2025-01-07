
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PTL.h"
#include "MB.h"
#include "EF.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,d;

  int numA,numB;
  int numatom;

  double length,f,crdA[3],crdB[3];

  double pi;
  
  char *line;
  size_t len=0;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *crdfilename,*parmfilename;

  FILE *crdfile,*parmfile,*outputfile,*logfile;
  
  char *progname;
  int opt_idx=1;
  
  pi=acos(-1.0);
  
  struct option long_opt[] = {
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
  
  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  crdfilename = *argv;
  parmfilename = *++argv;
  numA = atoi(*++argv);
  numB = atoi(*++argv);

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  if ( numA < 0 || numA > numatom  ) printf("error,there is no %d-th atom.\n",numA);
  if ( numB < 0 || numB > numatom  ) printf("error,there is no %d-th atom.\n",numB);

  crdfile=efopen(crdfilename,"r");
  getline(&line,&len,crdfile);
  fscanf(crdfile,"%d",&d);
  for (i=0;i<numatom;++i) {
    if (i==numA-1) for (j=0;j<3;++j) fscanf(crdfile,"%lf",&crdA[j]);
    else if (i==numB-1) for (j=0;j<3;++j) fscanf(crdfile,"%lf",&crdB[j]);
    else for (j=0;j<3;++j) fscanf(crdfile,"%lf",&f);
  }
  fclose(crdfile);

  length=0.0;
  for (i=0;i<3;++i) {
    length+=(crdA[i]-crdB[i])*(crdA[i]-crdB[i]);
  }
  length=sqrt(length);

  printf("%d-th %s to %d-th %s is %8.3lf\n",numA,AP.IGRAPH[numA-1],numB,AP.IGRAPH[numB-1],length);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] crdfilename parmfilename numA numB \n",progname);
}


