
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PT.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  int i,k;
  
  char *inputfilename1;
  FILE *inputfile1;
  
  if (argc < 2) {
    printf("USAGE: %s inputfilename(parmtop) \n",argv[0]);
    exit(1);
  }
  inputfilename1 = *++argv;
  inputfile1=efopen(inputfilename1,"r");
  readParmtop(inputfile1);
  fclose(inputfile1);

  for (k=0;k<AP.NPHIH;++k) {
    for (i=0;i<4;++i) {
      printf("%5d ",abs(AP.PH[k][i])/3+1);
    }
    printf("%5d ",abs(AP.PH[k][4]));
    printf("\n ");
  }

  for (k=0;k<AP.NPHIA;++k) {
    for (i=0;i<4;++i) {
      printf("%5d ",abs(AP.PA[k][i])/3+1);
    }
    printf("%5d ",abs(AP.PA[k][4]));
    printf("\n ");
  }
  
  return 0;
}

