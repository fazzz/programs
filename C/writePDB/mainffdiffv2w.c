#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PTLb.h"
#include "FFLc.h"

#include "EF.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,m;
  int numatomv,numatomw;

  int d;
  
  char *line;
  size_t len=0;

  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crdv,*crdw;

  struct potential ev,ew;
  struct force fv,fw;
  struct AmberParmL apv,apw;

  int c;

  char *crdvfilename,*crdwfilename,*parmfilenamev,*parmfilenamew,*progname;
  FILE *crdvfile,*crdwfile,*parmfilev,*parmfilew;

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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  crdvfilename = *argv;
  crdwfilename = *++argv;
  parmfilenamev = *++argv;
  parmfilenamew = *++argv;

  parmfilev=efopen(parmfilenamev,"r");
  readParmtopLb(parmfilev,&apv);
  fclose(parmfilev);
  numatomv=apv.NATOM;

  parmfilew=efopen(parmfilenamew,"r");
  readParmtopLb(parmfilew,&apw);
  fclose(parmfilew);
  numatomw=apw.NATOM;

  crdv=(double *)gcemalloc(sizeof(double)*numatomv*3);
  crdw=(double *)gcemalloc(sizeof(double)*numatomw*3);

  crdvfile=efopen(crdvfilename,"r");
  getline(&line,&len,crdvfile);
  fscanf(crdvfile,"%d",&d);
  for (i=0;i<numatomv;++i) for (j=0;j<3;++j) fscanf(crdvfile,"%lf",&crdv[i*3+j]);
  fclose(crdvfile);

  crdwfile=efopen(crdwfilename,"r");
  getline(&line,&len,crdwfile);
  fscanf(crdwfile,"%d",&d);
  for (i=0;i<numatomw;++i) for (j=0;j<3;++j) fscanf(crdwfile,"%lf",&crdw[i*3+j]);
  fclose(crdwfile);

  ffLc_set_calcffandforce(&ev,&fv,apv);
  ffLc_set_calcffandforce(&ew,&fw,apw);

  ffLc_calcffandforce(crdv,numatomv,&ev,&fv,apv);
  ffLc_calcffandforce(crdw,numatomw,&ew,&fw,apw);

  printf("%8.4e \n",ev.p_t-ew.p_t);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [--Amber] [-h] crdvfilename crdwfilename  parmtopvfilename parmtopwfilename\n",progname);
}
