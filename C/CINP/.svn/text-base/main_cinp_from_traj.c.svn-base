
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "CINP.h"
#include "PT.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j,amberflag,flag;
  int numsnap,numinpcrd,numinterval,numatom;
  int *num;

  char *option,*opt;  
  char *inputfilename1,*inputfilename2,*outputfilenamebase,*ext;
  FILE *inputfile1,*inputfile2;
  
  if (argc < 6) {
    printf("USAGE: ./%s [at] numsnap numinterval inputfilename1(trj) inputfilename2(parm)  outputfilenamebase(inp) \n",argv[0]);
    exit(1);
  }

  flag=(*++argv)[0];
  if (flag != 'a' && flag != 't') {
    printf("flag error: must be   or  ");
    exit(1);
  }
  if (flag=='a')
    amberflag=ON;
  else
    amberflag=OFF;
  numsnap=atoi(*++argv);
  numinterval=atoi(*++argv);
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  outputfilenamebase = *++argv;
  
  numinpcrd=(int)numsnap/numinterval;
  num=(int *)gcemalloc(sizeof(int)*numinpcrd);
  for (i=0;i<numinpcrd;++i)
    num[i]=i*numinterval;
  
  inputfile2=efopen(inputfilename2,"r");
  readParmtop(inputfile2);
  fclose(inputfile2);
  numatom=AP.NATOM;
  
  inputfile1=efopen(inputfilename1,"r");  
  Create_inpcrd_from_trj(outputfilenamebase,inputfile1,numsnap,numinpcrd,numatom,num,amberflag);
  fclose(inputfile1);  
  
  return 0;
}


