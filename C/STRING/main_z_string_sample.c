
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "CSI.h"
#include "EF.h"
#include "IO.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int flag;

  int numpoint,numiteration;

  double *p_stat/*,*p_evoluted,*force*/;

  double dt=0.01;

  char *line;
  size_t len=0;

  char *progname;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;

  progname=argv[0];

  while((c=getopt(argc,argv,"ht:"))!=-1) {
    switch(c) {
    case 't':
      dt=atof(optarg);
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  numiteration=atoi(*argv);
  numpoint=atoi(*++argv);
  inputfilename=*++argv;
  outputfilename=*++argv;

  //  p_stat = (double *)gcemalloc(sizeof(double)*numpoint*2);
  //  p_evoluted=(double *)gcemalloc(sizeof(double)*numpoint*2);
  //  force=(double *)gcemalloc(sizeof(double)*numpoint*2);

  p_stat = (double *)gcemalloc(sizeof(double)*numpoint*2);
  inputfile=efopen(inputfilename,"r");
  for (i=0;i<numpoint;++i) {
    fscanf(inputfile,"%lf",&p_stat[i*2]);
    fscanf(inputfile,"%lf",&p_stat[i*2+1]);
  }
  fclose(inputfile);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numiteration;++i) {
    z_string_sample(p_stat/*,p_evoluted,force*/,numpoint,dt);
    for (j=0;j<numpoint;++j)
      fprintf(outputfile,"%10.4lf %10.4lf\n",p_stat[j*2],p_stat[j*2+1]);
    fprintf(outputfile,"\n");
  }
  fclose(outputfile);

}

int USAGE(char *progname) {
  printf("USAGE:%s \n",progname);
  printf("-t dt    \n");
  printf("-h help  \n");
  printf("numpoint inputfilename outputfilename  \n");
}

 
