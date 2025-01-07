
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <netcdf.h>

#include "PDB.h"
#include "EF.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i;

  int numatom;
  PDBF PDBin,PDBout;

  char *progname;
  char *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"h"))!=-1) {
    switch(c) {
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  argc-=optind;
  argv+=optind;

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  outputfilename  = *++argv;

  inputfile=efopen(inputfilename,"r");
  readPDBatomnum(inputfile,&numatom);
  PDBin.numatom=numatom;
  PDBin.PDBa=gcemalloc(sizeof(PDBA)*numatom);
  inputfile=efopen(inputfilename,"r");
  readPDB(inputfile,PDBin,numatom);

  PDBout.numatom=numatom+2;
  PDBout.PDBa=gcemalloc(sizeof(PDBA)*(numatom+2));

  addCAP(PDBin,PDBout);

  outputfile=efopen(outputfilename,"w");
  writPDB(outputfile,PDBout);
  fclose(outputfile);

  return 0;
}

int USAGE(char *progname){
  printf("%s inputfilename outputfilename\n",progname);
}
