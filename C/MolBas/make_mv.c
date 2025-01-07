#include <stdio.h>
#include <stdlib.h>

#include "ParmTop.h"

int scandtraj(FILE *inputfile, double *coord, int numdihed);

int main(int argc, char *argv[]) {
  int i,j,k;
  int numdihed,numstep;
  double deltat,pi;
  double *dtraj,*Ftraj_R,*Ftraj_I,*PowerSpectral;
  char *inputfilename,*inputfilename2,*outputfilename;
  FILE *inputfile,*inputfile2, *outputfile;

  if (argc < 3) {
    printf("USAGE: ./make_mv.exe  inputfilename(traj) inputfilename2(parmtop)  outputfilename(pdb)\n");
    printf("cond: numstep numdihed deltat\n");
    exit(1);
  }
  inputfilename =  *++argv;
  inputfilename2 = *++argv;
  outputfilename = *++argv;

  readParmtop(inputfile);

  printf("MODEL\n");
  for(i=0;i<numlaststep;++i) {

    scandtraj(inputfile,coord,numdihed);
    for (natom=0;natom<AP.NATOM;++natom) {
      printf("ATOM%7d  %-2s           %10.6f%8.3f%8.3f\n",natom+1,AP.IGRAH[natom],coord[natom*3],coord[natom*3+1],coord[natom*3+2]);
    }
    printf("ENDMOD\n");
    printf("MODEL\n");
  }
}
