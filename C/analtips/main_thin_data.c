
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "IO.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j;
  int numinterval;
  
  char *inputfilename1,*outputfilename;
  FILE *inputfile1,*outputfile;
  
  if (argc < 4 ) {
    printf("USAGE: %s numinterval inputfilename1(data) outputfilename(thin_data)\n",argv[0]);
    exit(1);
  }
  numinterval=atoi(*++argv);
  inputfilename1 = *++argv;
  outputfilename = *++argv;
  
  inputfile1=efopen(inputfilename1,"r");
  outputfile=efopen(outputfilename,"w");
  io_thin_down_data(inputfile1,outputfile,numinterval);

  return 0;
}

