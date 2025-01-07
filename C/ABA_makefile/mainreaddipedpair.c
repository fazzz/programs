
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PT.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j;
  int flag=0;
  int *num;
  int **atomdihedpairs;
  
  char *inputfilename1,*inputfilename2,*outputfilename;
  FILE *inputfile1,*inputfile2, *outputfile;
  
  if (argc < 3) {
    printf("USAGE: %s flag(p or o or k) inputfilename1(pt) \n",argv[0]);
    exit(1);
  }
  flag=(*++argv)[0];
  if (flag != 'p' && flag != 'o' && flag != 'k') {
    printf("flag error: must be p or o or k\n");
    exit(1);
  }
  inputfilename1 = *++argv;
  outputfilename = *++argv;
  
  inputfile1=efopen(inputfilename1,"r");
  readParmtop(inputfile1);
  fclose(inputfile1);

  atomdihedpairs[0]=(int *)gcemalloc(sizeof(int)*4); // PHI,PSI
  atomdihedpairs[1]=(int *)gcemalloc(sizeof(int)*4); // OMEGA
  atomdihedpairs[2]=(int *)gcemalloc(sizeof(int)*4); // KI
  atomdihedpairs[3]=(int *)gcemalloc(sizeof(int)*8); // ACE
  atomdihedpairs[4]=(int *)gcemalloc(sizeof(int)*8); // NME
  num=(int *)gcemalloc(sizeof(int)*5);
  
  readdihedpairs(atomdihedpairs,num);

  for (i=0;i<num[3];++i) {
    printf("%4d-%4d-%4d-%4d\n",atomdihedpairs[3][i*4+0]+1,atomdihedpairs[3][i*4+1]+1,atomdihedpairs[3][i*4+2]+1,atomdihedpairs[3][i*4+3]+1);
  }
  for (i=0;i<num[0];++i) {
    printf("%4d-%4d-%4d-%4d\n",atomdihedpairs[0][i*4+0]+1,atomdihedpairs[0][i*4+1]+1,atomdihedpairs[0][i*4+2]+1,atomdihedpairs[0][i*4+3]+1);
  }
  for (i=0;i<num[1];++i) {
    printf("%4d-%4d-%4d-%4d\n",atomdihedpairs[1][i*4+0]+1,atomdihedpairs[1][i*4+1]+1,atomdihedpairs[1][i*4+2]+1,atomdihedpairs[1][i*4+3]+1);
  }
  for (i=0;i<num[2];++i) {
    printf("%4d-%4d-%4d-%4d\n",atomdihedpairs[2][i*4+0]+1,atomdihedpairs[2][i*4+1]+1,atomdihedpairs[2][i*4+2]+1,atomdihedpairs[2][i*4+3]+1);
  }
  for (i=0;i<num[4];++i) {
    printf("%4d-%4d-%4d-%4d\n",atomdihedpairs[4][i*4+0]+1,atomdihedpairs[4][i*4+1]+1,atomdihedpairs[4][i*4+2]+1,atomdihedpairs[4][i*4+3]+1);
  }
  
  return 0;
}

