
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PT.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  
  char *inputfilename1;
  FILE *inputfile1;
  
  if (argc < 2) {
    printf("USAGE: readnumatom inputfilename(parmtop) \n");
    exit(1);
  }
  inputfilename1 = *++argv;
  inputfile1=efopen(inputfilename1,"r");
  readParmtop(inputfile1);
  fclose(inputfile1);
  printf("%d\n",AP.NATOM);
  
  return 0;
}

