
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "IO.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j;
  int numinterval;
  
  char *inputfilename1,*inputfilename2,*outputfilename;
  FILE *inputfile1,*inputfile2, *outputfile;
  
  if (argc < 4) {
    printf("USAGE: ./thid_down.exe inputfilename1(data) inputfilename2(cond) outputfilename\n");
    exit(1);
  }
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  outputfilename = *++argv;
  
  inputfile2=efopen(inputfilename2,"r");
  fscanf(inputfile2,"%d",&numinterval);
  fclose(inputfile2);
  
  inputfile1=efopen(inputfilename1,"r");
  outputfile=efopen(outputfilename,"w");
  io_thin_down_data(inputfile1,outputfile,numinterval);
  fclose(inputfile1);
  fclose(outputfile);
  
  return 0;
}

