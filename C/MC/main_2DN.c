
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MB.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

int main(int argc, char *argv[]) {
  int i,j,k;
  char *line;
  size_t len=0;

  char *outputfilenam;
  FILE *outputfile;

  if (argc < 8) {
    printf("USAGE: ./%s \n",argv[0]);
    exit(1);
  }
  outputfilename = *++argv;

  logfile=efopen("log.txt","w");
  fprintf(logfile,"%d\n",numdihed);
  for (i=0;i<numdihed;++i) {
    for (j=0;j<4;++j) {
      fprintf(logfile,"%d ",atom_dihed_pair[i*4+j]);
    }
    for (j=0;j<4;++j) {
      fprintf(logfile,"%s ",AP.IGRAPH[atom_dihed_pair[i*4+j]-1]);
    }
    fprintf(logfile,"\n");
  }
  fclose(logfile);

  return 0;
}

