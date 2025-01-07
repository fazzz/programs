
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "Gaussian.h"
#include "EF.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,n;

  double x[2],minx=-1.0, maxx=1.0, dx=0.01, miny=-1.0, maxy=1.0, dy=0.01;
  double *nyu, **Sigma, pi;
  double *prob;

  char *outputfilename;
  FILE *outputfile;

  char *line;
  size_t len=0;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
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
      USAGE(progname);  exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < 1) {
    USAGE(progname);
    exit(1);
  }
  outputfilename = *argv;

  nyu=(double *)gcemalloc(sizeof(double)*2);
  nyu[0]=0.0;
  nyu[1]=0.0;

  Sigma=(double **)gcemalloc(sizeof(double *)*2);
  for (i=0;i<2;++i) Sigma[i]=(double *)gcemalloc(sizeof(double )*2);
  Sigma[0][0]=1.0;
  Sigma[0][1]=0.0;
  Sigma[1][0]=0.0;
  Sigma[1][1]=1.0;

  prob=Create_twoD_GaussianMap(minx, maxx, dx,
			       miny, maxy, dy,
			       nyu,  Sigma, pi);

  outputfile=efopen(outputfilename,"w");
  n=0;
  x[0]=minx;
  for (i=0;x[0]<maxx;++i) {
    x[0]=minx+dx*i;
    x[1]=miny;
    for (j=0;x[1]<maxy;++j){
      x[1]=miny+dy*j;
      fprintf(outputfile,"%8.3lf %8.3lf %8.3lf\n",x[0],x[1],prob[n]);
      ++n;
    }
  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename \n",progname);
}


