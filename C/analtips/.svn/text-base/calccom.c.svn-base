
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "MB.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"
#include "LA.h"

int main(int argc, char *argv[]) {
  int i,j,numstep,numatom;
  double *coord,*mass;
  double *com;
  
  char *inputfilename,*inputfilename2,*inputfilename3,*outputfilename;
  FILE *inputfile,*inputfile2,*inputfile3,*outputfile;
  
  if (argc < 5) {
    printf("USAGE: ./com.exe inputfilename(trj) inputfilename2(cond) inputfilename3(parmtop) outputfilename\n");
    exit(1);
  }
  inputfilename =  *++argv;
  inputfilename2 = *++argv;
  inputfilename3 = *++argv;
  outputfilename = *++argv;
  
  inputfile2=efopen(inputfilename2,"r");
  fscanf(inputfile2,"%d",&numstep);
  fclose(inputfile2);

  inputfile3=efopen(inputfilename3,"r");
  readParmtop(inputfile3);
  fclose(inputfile3);
  numatom=AP.NATOM;
  coord  = (double *)emalloc(sizeof(double)*numatom*3);
  mass  = (double *)emalloc(sizeof(double)*numatom);
  com  = (double *)emalloc(sizeof(double)*3);
  for (i=0;i<numatom;++i)
    mass[i]=AP.AMASS[i];
  
  inputfile=efopen(inputfilename,"r");
  outputfile=efopen(outputfilename,"w");
  for (i = 0; i < numstep; ++i) {
    io_scanconf(inputfile,numatom,coord,'x');
    MB_setcom(coord,mass,numatom,com);
    for (j = 0; j < 3; ++j)
      fprintf(outputfile,"%e ",com[j]);
    fprintf(outputfile,"\n");
  }
  
  fclose(inputfile);
  fclose(outputfile);

  free(coord);
  free(mass);
  free(com);
  
  return 0;
}

