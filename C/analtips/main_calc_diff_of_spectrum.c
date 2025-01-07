
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
  double f,dummy1;
  double freq,period,spect1,spect2;

  int sprow=3;
  double ifreq=0.0,ffreq=6000.0;

  int numstep,num;
  int numini=1,numfin,numrow=1;

  double criteria=0.1;

  double minperiod,maxfreq;
  int maxindex;

  char *inputfilename,*inputfilename2,*outputfilename,*progname;
  FILE *inputfile,*inputfile2,*outputfile;

  char *line,*dummy;
  size_t len=0;

  int c2;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c2=getopt(argc,argv,"hi:f:r:c:s:"))!=-1) {
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
    case 's':
      sprow=atoi(optarg);
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
  inputfilename2  = *++argv;
  outputfilename  = *++argv;

  inputfile=efopen(inputfilename,"r");
  inputfile2=efopen(inputfilename2,"r");
  outputfile=efopen(outputfilename,"w");
  for (i=0;;++i) {
    fscanf(inputfile,"%lf %lf",&freq,&period);
    fscanf(inputfile2,"%lf %lf",&freq,&period);
    fprintf(outputfile,"%lf %lf ",freq,period);
    if (freq < ifreq || freq > ffreq) break;
    for (j=2;j<numrow;++j)  {
      fscanf(inputfile,"%lf",&spect1);
      fscanf(inputfile2,"%lf",&spect2);
      fprintf(outputfile,"%10.8lf ",spect2-spect1);
    }
    fprintf(outputfile,"\n");
  }
  fclose(inputfile);
  fclose(inputfile2);
  fclose(outputfile);

  return 0;
}

void USAGE(char *progname) {
  printf("-i -- numini\n");
  printf("-f -- numfinp\n");
  printf("-r -- numrow\n");
  printf("-h -- help\n");
  printf("USAGE:%s [-i] [-f] [-j] [-r] [-h] inputfilename inputfilename2 outputfile\n", progname);
}
