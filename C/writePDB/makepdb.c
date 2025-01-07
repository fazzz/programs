#define _GNU_SOURCE  
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PDB.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

int io_inputtrj_Amberform2(FILE *inputfile,double *trj);

int main(int argc, char *argv[]) {
  int i,j,k,numatom;
  double *coord;
  PDBF PDB;
  
  char *inputfilename1,*inputfilename2,*outputfilename;
  FILE *inputfile1,*inputfile2, *outputfile;
  
  if (argc < 4) {
    printf("USAGE: makepdb inputfilename1(crd) inputfilename2(parm) outputfilename(pdb)\n");
    exit(1);
  }
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  outputfilename = *++argv;
  
  inputfile2=efopen(inputfilename2,"r");
  readParmtop(inputfile2);
  fclose(inputfile2);
  coord=(double *)emalloc(sizeof(double)*AP.NATOM*3);
  PDB.numatom=AP.NATOM;
  PDB.PDBa=(PDBA *)emalloc(sizeof(PDBA)*AP.NATOM);
  readPDBdatafromParmtop(PDB);
  
  inputfile1=efopen(inputfilename1,"r");
  io_inputtrj_Amberform2(inputfile1,coord);
  for (j=0;j<AP.NATOM;++j)
    for (k=0;k < 3; ++k)
      PDB.PDBa[j].coord[k]=coord[j*3+k];
  fclose(inputfile1);
  
  outputfile=efopen(outputfilename,"w");
  writPDB(outputfile,PDB);  
  fclose(outputfile);
  
  return 0;
}

int io_inputtrj_Amberform2(FILE *inputfile,double *trj) {
  int i,j;
  int numatom;
  char *line;
  size_t len=0;

  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d\n",&numatom);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      fscanf(inputfile,"%lf",&trj[i*3+j]);
    }
  }

  return 0;
}
