
#include <stdio.h>
#include <stdlib.h>

#include "PDB.h"
#include "EF.h"
#include "IO.h"

int main(int argc, char *argv[]) {
  int i,j,k,numatom,numstep;
  double *coord;
  PDBF PDB;
  char name[100][4];

  char *inputfilename,*inputfilename2,*outputfilename;
  FILE *inputfile,*inputfile2,*outputfile;
  
  if (argc < 4) {
    printf("USAGE: ./makemovie.exe inputfilename(trj) inputfilename2(cond) outputfilename\n");
    exit(1);
  }
  inputfilename =  *++argv;
  inputfilename2 = *++argv;
  outputfilename = *++argv;

  inputfile2=efopen(inputfilename2,"r");
  fscanf(inputfile2,"%d",&numstep);
  fscanf(inputfile2,"%d",&numatom);
  for (i=0;i<numatom;++i)
    fscanf(inputfile2,"%4s",&name[i]);
  fclose(inputfile2);

  coord=(double *)emalloc(sizeof(double)*numatom*3);
  PDB.numatom=numatom;
  PDB.PDBa=(PDBA *)emalloc(sizeof(PDBA)*numatom);

  for (i=0;i<numatom;++i) {
    PDB.PDBa[i].HETEROflag=0;
    PDB.PDBa[i].serial=i+1;
    PDB.PDBa[i].name[0]=' ';
    PDB.PDBa[i].name[1]=name[i][0];
    PDB.PDBa[i].name[2]=' ';
    PDB.PDBa[i].name[3]=' ';
    PDB.PDBa[i].altLOC=' ';
    PDB.PDBa[i].resname[0]=' ';
    PDB.PDBa[i].resname[1]=' ';
    PDB.PDBa[i].resname[2]=' ';
    PDB.PDBa[i].ChainID=' ';
    PDB.PDBa[i].resSeq=0;
    PDB.PDBa[i].iCode=' ';
    PDB.PDBa[i].occupancy=0.0;
    PDB.PDBa[i].tempfact=0.0;
  }

  inputfile=efopen(inputfilename,"r");
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    io_scanconf(inputfile,numatom,coord,'x');
    for (j=0;j<numatom;++j) {
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

