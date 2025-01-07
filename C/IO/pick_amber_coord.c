#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PT.h"
#include "IO.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j;
  int numatom,numstep;
  double *trj;
  
  char *inputfilename1,*inputfilename2,*inputfilename3,*outputfilename;
  FILE *inputfile1,*inputfile2,*inputfile3, *outputfile;
  
  if (argc < 5) {
    printf("USAGE: ./pick_amber_coord.exe inputfilename1(trj) inputfilename2(parmtop) inputfilename3(cond) outputfilename\n");
    exit(1);
  }
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  inputfilename3 = *++argv;
  outputfilename = *++argv;

  inputfile3=efopen(inputfilename3,"r");
  fscanf(inputfile3,"%d",&numstep);
  fclose(inputfile3);
  
  inputfile2=efopen(inputfilename2,"r");
  readParmtop(inputfile2);
  numatom=AP.NATOM;
  fclose(inputfile2);

  trj=(double *)ecalloc(sizeof(double),numatom*3);
  
  inputfile1=efopen(inputfilename1,"r");
  io_dismisstrj_Amberform(inputfile1,numstep-1,numatom);
  io_scanconf(inputfile1,numatom,trj,'x');
  fclose(inputfile1);
  
  outputfile=efopen(outputfilename,"w");
  io_outputconf_Amberform(outputfile,numatom,trj);
  fclose(outputfile);
  
  return 0;
}

