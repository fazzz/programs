
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j;
  int flag=0;
  
  char *inputfilename1,*inputfilename2,*inputfilename3,*outputfilename;
  FILE *inputfile1,*inputfile2,*inputfile2, *outputfile;
  
  if (argc < ) {
    printf("USAGE: %s flag() inputfilename1() inputfilename2() outputfilename\n",argv[0]);
    exit(1);
  }
  flag=(*++argv)[0];
  if (flag != ' ' && flag != ' ') {
    printf("flag error: must be   or  ");
    exit(1);
  }
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  inputfilename3 = *++argv;
  outputfilename = *++argv;
  
  inputfile2=efopen(inputfilename2,"r");
  
  fclose(inputfile2);
  
  inputfile1=efopen(inputfilename1,"r");
  
  fclose(inputfile1);
  
  outputfile=efopen(outputfilename,"w");
  
  fclose(outputfile);
  
  return 0;
}

