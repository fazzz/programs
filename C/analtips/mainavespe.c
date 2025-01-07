
#include <stdio.h>
#include <stdlib.h>

#include "PDB.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

int main(int argc, char *argv[]) {
  int i,j,k,d;
  int num;
  double f;
  double x[10000],spec[10000];

  char *inputfilename[10],*outputfilename;
  FILE *inputfile,*outputfile;
  
  for (i=0;i<10;++i) {
    inputfilename[i] =  *++argv;
  }
  outputfilename = *++argv;

  for (i=0;i<10;++i) {
    inputfile=efopen(inputfilename[i],"r");
    j=0;
    while (j<10000) {
      d=fscanf(inputfile,"%lf",&x[j]);
      d=fscanf(inputfile,"%lf",&f);
      spec[j]+=f/10.0;
      ++j;
    }
    fclose(inputfile);
  }

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<10000;++i) {
    fprintf(outputfile,"%lf %lf\n",x[i],spec[i]);
  }
  fclose(outputfile);

  return 1;
}
