//#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>

#include "ParmTop.h"

int main(int argc, char *argv[]){
  int i;
  char *inputfilename;
  FILE *inputfile;

  //  struct AmberParm AP;

  inputfilename=*++argv;

  if ((inputfile=fopen(inputfilename,"r"))==NULL){
    printf("error!");exit(1);
  }

  readParmtop(inputfile/*,AP*/);

  for (i=0;i<AP.NATOM;++i){
    printf("%16.8lf\n",AP.AMASS[i]);
  }

  return 1;
}
