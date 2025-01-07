
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include "EF.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  double pi;

  double ang,ange_ref,FC,dang;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;

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
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;  argc-=optind;  argv+=optind;

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  ange_ref       = atof(*argv);
  FC             = atof(*++argv);
  inputfilename  = *++argv;
  outputfilename = *++argv;

  inputfile=efopen(inputfilename,"r");
  outputfile=efopen(outputfilename,"w");

  while ((c=fscanf(inputfile,"%lf",&ang))!=-1){
    dang=(ang-ange_ref)*pi/180.0;
    if ( dang > pi ) dang-=2.0*pi;
    else if ( dang < -1.0*pi ) dang+=2.0*pi;
    fprintf(outputfile,"%20.10lf\n",0.5*FC*dang*dang);
  }

  fclose(inputfile);
  fclose(outputfile);

}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] ange_ref FC inputfilename outputfilename\n",progname);
}
