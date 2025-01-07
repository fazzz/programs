#define _GNU_SOURCE  
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "EF.h"
#include "IO.h"

int main(int argc, char *argv[]) {
  int i,j,k;
  int numtotalstep,numinterval,numvar,numcolum;
  int flag=0;
  char *line;
  size_t len=0;
  double *data,*sum,*sumofsq,*var,f;

  char *inputfilename,*inputfilename2,*inputfilename3,*outputfilename;
  FILE *inputfile,*inputfile2, *outputfile;
  
  if (argc < 4) {
    printf("USAGE: var flag(c or x) inputfilename(data) inputfilename2(cond) outputfilename\n");
    exit(1);
  }
  flag=(*++argv)[0];
  if (flag != 'c' && flag != 'x') {
    printf("flag error: must be c  or x ");
    exit(1);
  }
  inputfilename =  *++argv;
  inputfilename2 = *++argv;
  inputfilename3 = *++argv;
  outputfilename = *++argv;
  
  inputfile2=efopen(inputfilename2,"r");
  fscanf(inputfile2,"%d",&numtotalstep);
  fscanf(inputfile2,"%d",&numinterval);
  fscanf(inputfile2,"%d",&numcolum);
  fclose(inputfile2);
  
  numvar=(int)(numtotalstep/numinterval);
  data=(double *)ecalloc(sizeof(double),numcolum);
  sum=(double *)ecalloc(sizeof(double),numvar*numcolum);
  var=(double *)ecalloc(sizeof(double),numvar*numcolum);
  inputfile=efopen(inputfilename,"r");
  j=0;
  if (flag='c')
    getline(&line,&len,inputfile);
  for (i=0;i<numtotalstep;++i) {
    if ((i+1)%numinterval==0)
      ++j;
    io_scan_data(inputfile,numcolum,data,flag);
    for (k=0;k<numcolum;++k) {
       sum[j*numcolum+k]+=data[k];
       sumofsq[j*numcolum+k]+=data[k]*data[k];
    }
  }
  fclose(inputfile);
  for (i=0;i<numvar;++i)
    for (j=0;j<numcolum;++j)
      var[i*numcolum+j]=sumofsq[i*numcolum+j]/numinterval-(sum[i*numcolum+j]*sum[i*numcolum+j])/numinterval;
  
  outputfile=efopen(outputfilename,"w");
  io_outputtimeseries_f(outputfile,numvar,numcolum,var);
  fclose(outputfile);
  
   return 0;
}

