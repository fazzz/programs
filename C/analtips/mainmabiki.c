
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "PT.h"
#include "IO.h"
#include "EF.h"

int usage(void);

int main(int argc, char *argv[]) {
  int i,j,flag,flag2,c;
  int numlength,numdata,numatom;
  char *line;
  size_t len=0;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *inputfilename1,*inputfilename2,*outputfilename;
  FILE *inputfile1,*inputfile2,*outputfile;

  while((c=getopt(argc,argv,"LT:S:"))!=-1) {
    switch(c) {
    case 'T':
      flag='T';
      inputfilename2=optarg;
      break;
    case 'S':
      flag='S';
      numdata=atoi(optarg);
      break;
    case 'L':
      flag2='L';
      break;
    default:
      usage();
      exit(1);
    }
  }

  argc -= optind;

  argv += optind;

  if (argc < 3 ) {
    printf("USAGE: %s numlength inputfilename1(data) outputfilename(chop_data)\n",argv[0]);
    exit(1);
  }
  numlength=atoi(*argv);
  inputfilename1 = *++argv;
  outputfilename = *++argv;
  
  inputfile1=efopen(inputfilename1,"r");
  outputfile=efopen(outputfilename,"w");
  d1 = 1;
  while ( (c=getc(outputfile)) != -1  )  {
    if (c >= '0' || c <= '9') {

    }
    else if (c == '.' ) {

    }
    else if (c == '\b' || c == '\t') {

    }
    else if (c == '\n' ) {

    }
    else {

    }

  }
  fclose(inputfile1);
  fclose(outputfile);

  return 0;
}

int usage(void) {
  printf("USAGE:chop_data \n");
  printf("-T -S -L \n");
  printf("numlength inputfilename1 outputfilename\n");
}
