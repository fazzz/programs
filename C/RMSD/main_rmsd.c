
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "RMSD.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j;
  int flag=0;
  
  char *inputfilename1,*inputfilename2,*inputfilename3,*outputfilename;
  FILE *inputfile1,*inputfile2,*inputfile2, *outputfile;
  
  if (argc < 3) {
    printf("USAGE: %s flag(i or a) inputfilename1(trj) inputfilename2(cond) outputfilename(rmsd_trj)\n",argv[0]);
    exit(1);
  }
  flag=(*++argv)[0];
  if (flag != 'i' && flag != 'a') {
    printf("flag error: must be i  or a ");
    exit(1);
  }
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  outputfilename = *++argv;
  
  inputfile2=efopen(inputfilename2,"r");
  
  fclose(inputfile2);
  
  if (flag=='i') {
    inputfile1=efopen(inputfilename1,"r");
    outputfile=efopen(outputfilename,"w");
    io_scanconf(inputfile1,numatom,coordA,'x');
    for (i=0;i<numstep;++i) {
      io_scanconf(inputfile1,numatom,coordB,'x');
      rmsd=rmsd_qcp(coordA,coordB,mass,numatom);
    }
    fprintf(outputfile,"%d %e\n",i,rmsd);
    fclose(inputfile1);
    fclose(outputfile);
  }
  else {

  }
  
  
  return 0;
}

