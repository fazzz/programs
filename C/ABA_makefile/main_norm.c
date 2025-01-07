
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>

#include "IO.h"
#include "EF.h"
#include "mymath.h"

#define ON 1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,n;
  double f;

  int numstep,num;
  int numini=1,numfin,numinirow,numfinrow,numrow;

  double *dat,*dat_norm;

  char *inputfilename,*outputfilename,*progname;
  FILE *inputfile,*outputfile;

  char *line,*dummy;
  size_t len=0;

  int c2;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c2=getopt(argc,argv,"hi:f:j:g:r:"))!=-1) {
    switch(c2) {
    case 'i':
      numini=atof(optarg)-1;
      break;
    case 'f':
      numfin=atoi(optarg)-1;
      break;
    case 'j':
      numinirow=atoi(optarg)-1;
      break;
    case 'g':
      numfinrow=atoi(optarg)-1;
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

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  outputfilename = *++argv;

  inputfile=efopen(inputfilename,"r");
  for (i=0;i<numini;++i)  getline(&line,&len,inputfile);

  numstep=(int)(numfin-numini+1);
  num=numfinrow-numinirow+1;

  dat=(double *)gcemalloc(sizeof(double)*numstep);
  dat_norm=(double *)gcemalloc(sizeof(double)*numstep);

  for (i=0;i<numstep;++i) 
    for (j=0;j<numrow;++j) 
      if (j>=numinirow && j<=numfinrow ) fscanf(inputfile,"%lf",&dat[i]);
      else fscanf(inputfile,"%lf",&f);

  dat_norm=ts_normalize_od(dat,numstep);
  fclose(inputfile);
    
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) fprintf(outputfile,"%e \n",dat_norm[i]);
  fclose(outputfile);

  return 0;
}

void USAGE(char *progname) {
  printf("-i -- numini\n");
  printf("-f -- numfinp\n");
  printf("-j -- numinirow\n");
  printf("-g -- numfinrow\n");
  printf("-r -- numrow\n");
  printf("-h -- help\n");
  printf("USAGE:%s [-i] [-j] [-g] [-r] [-h] inputfilename outputfile\n", progname);
}
