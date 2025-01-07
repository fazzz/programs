#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PTL.h"
#include "FFL.h"

#include "EF.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,m;
  int numatomv;

  int d;
  
  char *line;
  size_t len=0;

  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crdv,*crdw;

  struct potential ev;
  struct force fv;

  int c;

  char *crdvfilename,*parmfilenamev,*progname;
  FILE *crdvfile,*parmfilev;

  int opt_idx=1;

  struct option long_opt[] = {
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"h",long_opt,&opt_idx))!=-1) {
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
  crdvfilename = *argv;
  parmfilenamev = *++argv;

  parmfilev=efopen(parmfilenamev,"r");
  readParmtopL(parmfilev);
  fclose(parmfilev);
  numatomv=AP.NATOM;

  crdv=(double *)gcemalloc(sizeof(double)*numatomv*3);

  crdvfile=efopen(crdvfilename,"r");
  getline(&line,&len,crdvfile);
  fscanf(crdvfile,"%d",&d);
  for (i=0;i<numatomv;++i) for (j=0;j<3;++j) fscanf(crdvfile,"%lf",&crdv[i*3+j]);
  fclose(crdvfile);

  ffL_set_calcffandforce(&ev,&fv);
  
  ffL_calcffandforce(crdv,numatomv,&ev,&fv);

  printf("%8.4e \n",ev.p_t);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [--Amber] [-h] crdvfilename parmtopvfilename \n",progname);
}
