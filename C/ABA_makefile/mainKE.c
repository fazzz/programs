
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "EF.h"
#include "PT.h"

#define kbkcl 1.98723e-3
#define kbuap 1.98723e-3*4.18407*100.0

int main(int argc, char *argv[]) {
  int i,j;
  int numatom,numstep;
  int flag=0;
  double KE,*trj;

  char *line;
  size_t len=0;
  
  char *inputfilename1,*inputfilename2,*outputfilename;
  FILE *inputfile1,*parmtop,*outputfile;
  
  if (argc < 4) {
    printf("USAGE: %s flag(a or t) numstep inputfilename1(trj) inputfilename2(parmtop) \n",argv[0]);
    exit(1);
  }
  flag=(*++argv)[0];
  if (flag != 'a' && flag != 't') {
    printf("flag error: must be a  or t ");
    exit(1);
  }
  numstep=atoi(*++argv);
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  
  parmtop=efopen(inputfilename2,"r");
  readParmtop(parmtop);
  fclose(parmtop);

  numatom=AP.NATOM;
  inputfile1=efopen(inputfilename1,"r");
  trj=(double *)ecalloc(sizeof(double),numstep*numatom*3);
  if (flag=='a')
    getline(&line,&len,inputfile1);  
  io_scantraj_aw(inputfile1,numstep,numatom,trj);
  fclose(inputfile1);

  KE=0.0;
  for (i=0;i<numstep*numatom*3;++i)
    KE+=0.5*trj[i]*trj[i]/numstep;

  printf("%lf12.4 kcal/mol\n",KE);
  printf("%lf12.4 K\n",KE/((3*numatom)-6)/*-(AP.NBONA-AP.MBONA)-(AP.NTHETA-AP.MTHETA)-(AP.NPHIA-AP.MPHIA)*//kbkcl/2);
  
  return 0;
}

