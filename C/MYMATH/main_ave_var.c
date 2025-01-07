
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "EF.h"
#include "IO.h"

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,d,n;
  double f;
  int nrow=1,interval=1,nini=0;
  double *ave,*var;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *progname;
  char *inputfilename;
  char *outputfilename;

  FILE *inputfile;
  FILE *outputfile;

  while((c=getopt(argc,argv,"hn:i:j:"))!=-1) {
    switch(c) {
    case 'n':
      nrow=atoi(optarg);
      break;
    case 'i':
      interval=atoi(optarg);
      break;
    case 'j':
      nini=atoi(optarg);
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

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  outputfilename = *++argv;

  ave=(double *)gcemalloc(sizeof(double)*nrow);
  var=(double *)gcemalloc(sizeof(double)*nrow);

  for (i=0;i<nrow;++i) {
    ave[i]=0.0;
    var[i]=0.0;
  }

  inputfile=efopen(inputfilename,"r");
  n=0;
  i=0;
  while ( d!= -1 ) {
    if (i<nini)
      getline(&line,&len,inputfile);
    else {
      for (j=0;j<nrow;++j) {
	d=fscanf(inputfile,"%lf",&f);
	if (i%interval==0 ) {
	  //	  ave[j]=(n*ave[j]+f)/(n+1);
	  //	  var[j]=(n*var[j]+f*f)/(n+1);
	  ave[j]+=f;
	  var[j]+=f*f;
	  if (j==0) ++n;
	}
      }
    }
    ++i;
  }
  fclose(inputfile);

  for (i=0;i<nrow;++i) ave[i]/=n;
  for (i=0;i<nrow;++i) var[i]/=n;
  for (i=0;i<nrow;++i) var[i]-=ave[i]*ave[i];
  for (i=0;i<nrow;++i) var[i]=sqrt(var[i]);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<nrow;++i) fprintf(outputfile,"%e %e\n",ave[i],var[i]);
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s inputfilename outputfilename \n",progname);
}
