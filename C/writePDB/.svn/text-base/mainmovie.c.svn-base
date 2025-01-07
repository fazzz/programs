
#include <stdio.h>
#include <stdlib.h>

#include "PDB.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

int main(int argc, char *argv[]) {
  int i,j,k,numstep;
  double *coord;
  PDBF PDB;

  char *inputfilename,*inputfilename2,*inputfilename3,*outputfilename;
  FILE *inputfile,*inputfile2,*inputfile3,*outputfile;
  
  if (argc < 5) {
    printf("USAGE: ./makemovie.exe numstep inputfilename(trj) inputfilename3(parmtop) outputfilename\n");
    exit(1);
  }
  numstep=atoi(*++argv);
  inputfilename =  *++argv;
  inputfilename3 = *++argv;
  outputfilename = *++argv;

  inputfile3=efopen(inputfilename3,"r");
  readParmtop(inputfile3);
  fclose(inputfile3);
  coord=(double *)emalloc(sizeof(double)*AP.NATOM*3);
  PDB.numatom=AP.NATOM;
  PDB.PDBa=(PDBA *)emalloc(sizeof(PDBA)*AP.NATOM);
  readPDBdatafromParmtop(PDB);

  inputfile=efopen(inputfilename,"r");
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    io_scanconf(inputfile,AP.NATOM,coord,'x');
    for (j=0;j<AP.NATOM;++j) {
      for (k=0;k < 3; ++k)	{ 
	PDB.PDBa[j].coord[k]=coord[j*3+k];
      }
    }
    fprintf(outputfile,"MODEL\n");
    writPDB(outputfile,PDB);
    fprintf(outputfile,"ENDMOD\n");
  }
  fclose(inputfile);
  fclose(outputfile);

  return 1;
}
