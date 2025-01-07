
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PT.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  double s;
  
  char *inputfilename1,*outputfilename;
  FILE *inputfile1,*outputfile;
  
  if (argc < 3) {
    printf("USAGE: scale(double) inputfilename1(pt) outputfilename(ptnew)\n",argv[0]);
    exit(1);
  }
  s=atof(argv[1]);
  inputfilename1 = argv[2];
  outputfilename = argv[3];
  
  inputfile1=efopen(inputfilename1,"r");
  readParmtop(inputfile1);
  fclose(inputfile1);

  sd_woba(s);

  outputfile=efopen(outputfilename,"w");
  writeParmtop(outputfile);
  fclose(outputfile);

}

