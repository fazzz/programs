#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "PT.h"

void usage(char *progname);

#define ON 1
#define OFF 0

int main(int argc, char *argv[]) {
  int i,j;
  int outflag=OFF;
  int nres;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *parmfilename,*outputfilename,*progname;
  FILE *parmfile,*outputfile;

  while((c=getopt(argc,argv,"hf:"))!=-1) {
    switch(c) {
    case 'f':
      outflag=ON;
      outputfilename=optarg;
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
  parmfilename = *argv;

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);

  nres=AP.NRES;

  if (outflag==OFF) {
    for (i=0;i<nres;++i) {
      printf("%s ",AP.LABERES[i]);
    }
    printf("\n");
  }
  else {
    outputfile=efopen(outputfilename,"w");
    for (i=0;i<nres;++i) 
      fprintf(outputfile,"%s ",AP.LABERES[i]);
    fprintf(outputfile,"\n");
    fclose(outputfile);
  }
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-f outputfilename  ] outputfilemode \n");
  printf("[-h] help \n");
  printf("%s [-f outputfilename ]  [-h] parmfilename \n",progname);
}
