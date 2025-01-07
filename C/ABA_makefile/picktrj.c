
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "IO.h"
#include "PT.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j;
  int numstep,numatom,interval;
  double **trj;
  char *line;
  size_t len=0;
  
  char *inputfilename1,*inputfilename2,*outputfilename;
  FILE *inputfile1,*inputfile2,*outputfile;
  
  if (argc < 6) {
    printf("USAGE: %s  totalstep interval inputfilename1(crd)  inputfilename2(parm) outputfilename1(pc) \n",argv[0]);
    exit(1);
  }
  numstep=atof(*++argv);
  interval=atof(*++argv);
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  outputfilename = *++argv;
  
  inputfile2=efopen(inputfilename2,"r");
  readParmtop(inputfile2);
  fclose(inputfile2);
  numatom=AP.NATOM;

  trj=(double **)gcemalloc(sizeof(double *));
  trj[0]=(double *)gcemalloc(sizeof(double)*(numatom*3));
  
  inputfile1=efopen(inputfilename1,"r");
  outputfile=efopen(outputfilename,"w");

  for (i=0;i<numstep;++i) {
    io_scantrj(inputfile1,numatom,1,trj);

    if (i%interval==0)
      io_outtrj(outputfile,numatom,1,trj);
  }
  
  fclose(inputfile1);
  fclose(inputfile2);
  fclose(outputfile);

  return 0;
}

