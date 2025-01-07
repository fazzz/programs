#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"
#include "GOLMAA_MB_PROTEINS2008.h"
#include "GOLMAA_MB_PROTEINS2008_check.h"

#include "PTL.h"
#include "EF.h"
#include "NC.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k;

  double xi=-10,xf=10;
  double x1=-1.0,x2=1.0;
  double xd=0.001;

  double x;
  double e,e1,e2,*f;
  double d=1.0,de=1.0,d2;

  double k1=10.0,k2=10.0;

  double kai;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *outputfilename;
  FILE *outputfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"xi",1,NULL,'i'},
    {"xf",1,NULL,'f'},
    {"xd",1,NULL,'p'},
    {"k1",1,NULL,'1'},
    {"k2",1,NULL,'2'},
    {"de",1,NULL,'d'},
    {"d",1,NULL,'x'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"i:f:p:1:2:d:x:h",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'i':
      xi=atof(optarg);
      break;
    case 'f':
      xf=atof(optarg);
      break;
    case 'p':
      xd=atof(optarg);
      break;
    case '1':
      k1=atof(optarg);
      break;
    case '2':
      k2=atof(optarg);
      break;
    case 'd':
      de=atof(optarg);
      break;
    case 'x':
      d=atof(optarg);
      break;
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

  if (argc < 1) {
    USAGE(progname);
    exit(1);
  }
  outputfilename = *argv;

  d2=d*d;
  f=(double *)gcemalloc(sizeof(double)*3);

  outputfile=efopen(outputfilename,"w");
  for (x=xi;x<xf;x+=xd) {
    e=GOLMAA_MB_PROTEINS2008_ff_calcff_harmo(x,de,d2,k1,x1,k2,x2,&f);
    kai=GOLMAA_MB_PROTEINS2008_harmo_Kai(x,de,d,d2,k1,x1,k2,x2);

    e1=0.5*k1*(x-x1)*(x-x1);
    e2=0.5*k2*(x-x2)*(x-x2);

    fprintf(outputfile,"%e %e %e %e %e\n",x,e,e1,e2,kai);
  }
  fclose(outputfile);


  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfilename parmfilename outputfilename\n",progname);
}

