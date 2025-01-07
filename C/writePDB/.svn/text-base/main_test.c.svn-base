

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "EF.h"
#include "PDB.h"

int main(int argc, char *argv[]) {
  int numatom;
  char *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;

  PDBF PDBin,PDBout;
  
  if (argc < 3) {
    printf("USAGE: ./PDB_test.exe inputfilename(PDB)  outputfilename(PDB)\n");
    exit(1);
  }

  inputfilename =  *++argv;
  outputfilename = *++argv;
  
  inputfile=efopen(inputfilename,"r");
  //  readPDBatomnum(inputfile,numatom);
  // fclose(inputfile);
  inputfile=efopen(inputfilename,"r");
  readPDB(inputfile,PDBin,20);
  fclose(inputfile);
  
  copyPDBform(PDBin,PDBout);

  outputfile=efopen(outputfilename,"w");
  writePDB(outputfile,PDBout);
  fclose(outputfile);
  
  return 0;
}

