
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>

#include "IO.h"
#include "EF.h"

#define ON 1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,n;
  double *f;
  double freq,period;

  double ifreq=0.0,ffreq=6000.0;

  int numstep,num;
  int numini=1,numfin,numrow;

  double criteria=0.1;

  double minperiod,maxfreq;
  int maxindex;

  char *inputfilename,*outputfilename,*progname;
  FILE *inputfile,*outputfile;

  char *line,*dummy;
  size_t len=0;

  int c2;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c2=getopt(argc,argv,"hi:f:r:c:"))!=-1) {
    switch(c2) {
    case 'i':
      ifreq=atof(optarg);
      break;
    case 'f':
      ffreq=atof(optarg);
      break;
    case 'r':
      numrow=atoi(optarg);
      break;
    case 'c':
      criteria=atof(optarg);
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

  if (argc < 1) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  outputfilename = *++argv;

  inputfile=efopen(inputfilename,"r");
  for (i=0;i<numini;++i)  getline(&line,&len,inputfile);

  f=(double *)gcemalloc(sizeof(double)*(numrow-2));

  for (i=0;;++i) {
    fscanf(inputfile,"%lf",&freq);
    if (freq <= ifreq || freq > ffreq) break;
    fscanf(inputfile,"%lf",&period);    
    if (i==0) minperiod=period;
    for (j=2;j<numrow;++j) {
      fscanf(inputfile,"%lf",&f[j-2]);
      if ( minperiod > period && f[j-2] > criteria) {
	maxindex=j;
	maxfreq=freq;
	minperiod=period;
      }
    }
  }
  fclose(inputfile);

  outputfile=efopen(outputfilename,"w");
  fprintf(outputfile,"%d %e %e\n",maxindex+1,minperiod,maxfreq);
  fclose(outputfile);

  return 0;
}

void USAGE(char *progname) {
  printf("-i -- numini\n");
  printf("-f -- numfinp\n");
  printf("-r -- numrow\n");
  printf("-h -- help\n");
  printf("USAGE:%s [-i] [-f] [-j] [-r] [-h] inputfilename outputfile\n", progname);
}
