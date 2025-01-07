#include <stdio.h>
#include <stdlib.h>

#include "PDB.h"

int main(int argc, char *argv[]) {
  char *inputfilename,*outputfilename;

  if (argc < 2) {
    printf("USAGE: ./checkwrirePDB.exe inputPDB outputPDB\n");
    exit(1);
  }
  inputfilename =  *++argv;
  outputfilename = *++argv;

  readPDB(inputfilename);
  writePDB(outputfilename);

  return 1;
}
