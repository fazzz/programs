
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SSL.h"
#include "EF.h"
#include "IO.h"

int outputtimeseries(FILE *outputfile,int numstep,int numdimension,double *data,int flag);

int main(int argc, char *argv[]) {
  int i,j,k,flag;
  int nums,numv;
  double *data,*datanorm;

  char *inputfilename,*inputfilename2;
  char *outputfilename;
  FILE *inputfile,*inputfile2;
  FILE *outputfile;

  if (argc < 3) {
    printf("USAGE: ./sslnorm.exe flag (c or v) inputfilename(data) inputfilename2(cond) outputfilename\n");
    exit(1);
  }
  flag=(*++argv)[0];
  if (flag != 'c' && flag != 'v') {
    printf("flag error: must be c or v");
    exit(1);
  }

  inputfilename   = *++argv;
  inputfilename2  = *++argv;
  outputfilename  = *++argv;

  inputfile2=efopen(inputfilename2,"r");
  fscanf(inputfile2,"%d",&nums);
  fscanf(inputfile2,"%d",&numv);
  fclose(inputfile2);

  inputfile   =  efopen(inputfilename,"r");
  outputfile  = efopen(outputfilename,"w");

  data=(double *)ecalloc(numv*nums,sizeof(double));
  datanorm=(double *)ecalloc(numv*nums,sizeof(double));

  io_scantimeseries(inputfile,nums,numv,data,flag);
  ssl_normalize(data,nums,numv,datanorm);
  outputtimeseries(outputfile,nums,numv,datanorm,'c');

  fclose(inputfile);
  fclose(outputfile);
  free(data);
  free(datanorm);

  return 0;
}

int outputtimeseries(FILE *outputfile,int numstep,int numdimension,double *data,int flag){
  int i,j,d;
  
  for (i=0;i<numstep;++i) {
    if (flag=='c')
      fprintf(outputfile,"%d ",i);
    for (j=0;j<numdimension;++j)
      fprintf(outputfile,"%10.6lf",data[i*numdimension+j]);
    fprintf(outputfile,"\n");
  }
  fprintf(outputfile,"\n");
  
  return 0;
}
