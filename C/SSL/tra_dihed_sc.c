
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "efunc.h"

int main(int argc, char *argv[]) {
  int i,j;
  int nums,numv;
  double f,pi;

  char *inputfilename,*inputfilename2,*outputfilename;
  FILE *inputfile,*inputfile2, *outputfile;

  if (argc < 3) {
    printf("USAGE: ./tds.exe inputfilename(data) inputfilename2(cond) outputfilename\n");
    exit(1);
  }
  inputfilename =  *++argv;
  inputfilename2 = *++argv;
  outputfilename = *++argv;

  inputfile2=efopen(inputfilename2,"r");
  fscanf(inputfile2,"%d",&nums);
  fscanf(inputfile2,"%d",&numv);
  fclose(inputfile2);

  inputfile=efopen(inputfilename,"r");
  outputfile=efopen(outputfilename,"w");

  pi=acos(-1.0);
  for (i=0;i<nums;++i) {
    for (j=0;j<numv;++j) {
      fscanf(inputfile,"%lf",&f);
      f = (f + 180.0)/pi;
      fprintf(outputfile,"%e %e ",cos(f),sin(f));
    }
    fprintf(outputfile,"\n");
  }

  fclose(inputfile);
  fclose(outputfile);

}

