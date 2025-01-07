
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j;
  int flag=0;

  int num,numdata;
  double width_orign,width,sum;
  double *data;
  
  char *inputfilename1,*outputfilename;
  FILE *inputfile1,*outputfile;
  
  if (argc < 5) {
    printf("USAGE: %s width numdata inputfilename1(data) outputfilename\n",argv[0]);
    exit(1);
  }
  width=atoi(*++argv);
  numdata=atoi(*++argv);
  inputfilename1 = *++argv;
  outputfilename = *++argv;

  data=(double *)gcemalloc(sizeof(double)*numdata*2);
  
  inputfile1=efopen(inputfilename1,"r");
  outputfile=efopen(outputfilename,"w");
  io_scnasdcoldata(inputfile1,numdata,0,2,1,2,data);
  num=0;

  //  width_orign=data[(i+1)*2]-data[i*2];
  for (i=0;i<numdata;++i) {
    if (data[i*2] < data[0]+width*num) {
      sum+=data[i*2+1];
    }
    else {
      sum=sum/width;
      fprintf(outputfile,"%e %e \n",data[0]+(width*(num-0.5)),/*(int)*/sum);
      sum=data[i*2+1];
      num+=1;
    }
  }

  fclose(inputfile1);
  fclose(outputfile);
  
  return 0;
}

