
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "writeAmberInput.h"

int writAmberInput(FILE *inputfile,int  numatom, double *crd, char NAME[3]) {  
  int i;

  for (i=0;i<3;++i) fprintf(inputfile,"%c",NAME[i]);  fprintf(inputfile,"\n");

  fprintf(inputfile,"%5d\n",numatom);

  for (i=0;i<numatom*3;++i) {
    fprintf(inputfile,"%11.7lf ",crd[i]);

    if ((i+1)%6==0) fprintf(inputfile,"\n");
  }
  return 1;
}
