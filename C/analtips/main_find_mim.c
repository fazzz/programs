#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EF.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,n=0;

  int numbinx,numbiny;

  double **x,**y,**pmf;
  
  char *inputfilename,*outputfilename_base,outputfilename[1000];
  FILE *inputfile, *outputfile;

  char *line;
  size_t len=0;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;
  
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
  
  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  numbinx = atoi(*argv);
  numbiny = atoi(*++argv);
  inputfilename = *++argv;
  outputfilename_base = *++argv;

  x=(double **)gcemalloc(sizeof(double *)*numbinx);
  for (i=0;i<numbinx;++i) x[i]=(double *)gcemalloc(sizeof(double)*numbiny); 
  y=(double **)gcemalloc(sizeof(double *)*numbinx);
  for (i=0;i<numbinx;++i) y[i]=(double *)gcemalloc(sizeof(double)*numbiny); 
  pmf=(double **)gcemalloc(sizeof(double *)*numbinx);
  for (i=0;i<numbinx;++i) pmf[i]=(double *)gcemalloc(sizeof(double)*numbiny); 

  inputfile=efopen(inputfilename,"r");
  for (i=0;i<numbinx;++i) {
    for (j=0;j<numbiny;++j) {
      fscanf(inputfile,"%lf",&x[i][j]);
      fscanf(inputfile,"%lf",&y[i][j]);
      fscanf(inputfile,"%lf",&pmf[i][j]);
    }
  }
  fclose(inputfile);

  for (i=0;i<numbinx;++i) {
    for (j=0;j<numbiny;++j) {
      if ( i-1 > 0 && j-1 > 0 && i+1 < numbinx && j+1 < numbiny ) {
	if ( pmf[i][j] < pmf[i+1][j]
	     && pmf[i][j] < pmf[i][j+1]
	     && pmf[i][j] < pmf[i+1][j+1]
	     && pmf[i][j] < pmf[i-1][j]
	     && pmf[i][j] < pmf[i][j-1]
	     && pmf[i][j] < pmf[i-1][j-1]
	     && pmf[i][j] < pmf[i+1][j-1]
	     && pmf[i][j] < pmf[i-1][j+1]
	     && pmf[i][j] < 100.0) {
	  printf("%d %8.6lf %8.6lf %8.6lf\n",n+1,x[i][j],y[i][j],pmf[i][j]);
	  ++n;
	}
      }
    }
  }

  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] numbinx numbiny inputfilename \n",progname);
}


