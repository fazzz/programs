
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
  double f;

  int numstep,num;
  int numini=1,numfin,numspecrow,numrow;

  double max;
  int maxindex;

  char *inputfilename,*outputfilename,*progname;
  FILE *inputfile,*outputfile;

  char *line,*dummy;
  size_t len=0;

  int c2;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c2=getopt(argc,argv,"hi:f:j:r:"))!=-1) {
    switch(c2) {
    case 'i':
      numini=atof(optarg)-1;
      break;
    case 'f':
      numfin=atoi(optarg)-1;
      break;
    case 'j':
      numspecrow=atoi(optarg)-1;
      break;
    case 'r':
      numrow=atoi(optarg);
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

  numstep=(int)(numfin-numini+1);

  inputfile=efopen(inputfilename,"r");
  for (i=0;i<numini;++i)  getline(&line,&len,inputfile);

  for (i=0;i<numstep;++i) {
    for (j=0;j<numrow;++j) {
      if (j==numspecrow) {
	fscanf(inputfile,"%lf",&f);
	if (i==0) {
	  max=f;
	  maxindex=i;
	}
	if ( max < f) {
	  max=f;
	  maxindex=i;
	}
      }
      else fscanf(inputfile,"%lf",&f);
    }
  }
  fclose(inputfile);

  inputfile=efopen(inputfilename,"r");
  for (i=0;i<maxindex+1;++i)  getline(&line,&len,inputfile);
  fclose(inputfile);
  
  printf("%s\n",line);

  return 0;
}

void USAGE(char *progname) {
  printf("-i -- numini\n");
  printf("-f -- numfinp\n");
  printf("-j -- numspecrow\n");
  printf("-r -- numrow\n");
  printf("-h -- help\n");
  printf("USAGE:%s [-f]  [-j] [-r] [-h] inputfilename outputfile\n", progname);
}
