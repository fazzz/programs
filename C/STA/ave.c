#define _GNU_SOURCE  
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "EF.h"
#include "IO.h"

int main(int argc, char *argv[]) {
  int i,j,k,d;
  int numtotalstep,numinterval,numave,numcolum;
  int flag=0;
  double *data,*sum,*ave,f;
  
  char *inputfilename,*inputfilename2,*outputfilename;
  FILE *inputfile,*inputfile2, *outputfile;
  
  if (argc < 4) {
    printf("USAGE: ave flag(c or x) inputfilename(data) inputfilename2(cond) outputfilename(ave)\n");
    printf("cond: numtotalstep numinterval numcolum\n");
    exit(1);
  }
  flag=(*++argv)[0];
  if (flag != 'c' && flag != 'x') {
    printf("flag error: must be c   or x ");
    exit(1);
  }
  inputfilename =  *++argv;
  inputfilename2 = *++argv;
  outputfilename = *++argv;
  
  inputfile2=efopen(inputfilename2,"r");
  fscanf(inputfile2,"%d",&numtotalstep);
  fscanf(inputfile2,"%d",&numinterval);
  fscanf(inputfile2,"%d",&numcolum);
  fclose(inputfile2);

  numave=(int)(numtotalstep/numinterval);
  data=(double *)egmalloc(sizeof(double),numcolum);
  sum=(double *)egmalloc(sizeof(double),numave*numcolum);
  ave=(double *)egmalloc(sizeof(double),numave*numcolum);
  inputfile=efopen(inputfilename,"r");
  j=0;
  if (flag='c')
    getline(&line,&len,inputfile);
  for (i=0;i<numtotalstep;++i) {
    if ((i+1)%numinterval==0)
      ++j;
    io_scan_data(inputfile,numcolum,data,flag);
    for (k=0;k<numcolum;++k)
      sum[j*numcolum+k]+=data[k];
  }
  fclose(inputfile);
  for (i=0;i<numave;++i)
    for (j=0;j<numcolum;++j)
    ave[i*numcolum+j]=sum[i*numcolum+j]/numinterval;
  
  outputfile=efopen(outputfilename,"w");
  io_outputtimeseries_f(outputfile,numave,numcolum,ave);
  fclose(outputfile);
  
  return 0;
}
