
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PT.h"
#include "IO.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j,numatom;
  double KE;
  double *velo;
  
  char *inputfilename1,*inputfilename2;
  FILE *inputfile1,*inputfile2;
  
  if (argc < 3) {
    printf("USAGE: %s inputfilename1(vel) inputfilename2(parmtop)  \n",argv[0]);
    exit(1);
  }
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  
  inputfile2=efopen(inputfilename2,"r");
  readParmtop(inputfile2);  
  fclose(inputfile2);
  numatom=AP.NATOM;

  velo=(double *)ecalloc(sizeof(double),numatom*3);
  inputfile1=efopen(inputfilename1,"r");
  io_scanconf(inputfile1,numatom,velo,'x');
  fclose(inputfile1);

  KE=0.0;
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j)
      KE+=0.5*velo[i*3+j]*20.455*velo[i*3+j]*20.455*AP.AMASS[i]/(4.18407*100.0);

  printf("%e\n",KE);
  for (i=0;i<numatom;++i)
    printf("%e %e %e \n",velo[i*3],velo[i*3+1],velo[i*3+2]);

  free(velo);
  
  return 0;
}

