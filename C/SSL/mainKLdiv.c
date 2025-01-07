
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SSL.h"
#include "EF.h"
#include "IO.h"

int outputKLdiv(FILE *outputfile,double *KLdiv,int numv);
int outputKLdivtopten(FILE *outputfile,int topten[10], double vtopten[10], int num);

int main(int argc, char *argv[]) {
  int i,j,k,num;
  int numv,topten[10];
  double *Lambda,*Sigma,*Lambdaini,*Sigmaini,*KLdiv,vtopten[10];

  char *inputfilename,*inputfilename2,*inputfilename3;
  char *outputfilename,*outputfilename2;
  FILE *inputfile,*inputfile2,*inputfile3;
  FILE *outputfile,*outputfile2;

  if (argc <5) {
    printf("USAGE: ./KLdiv.exe inputfilename1(data1) inputfilename2(data2) inputfilename3(cond) outputfilename\n");
    exit(1);
  }

  inputfilename   = *++argv;
  inputfilename2  = *++argv;
  inputfilename3  = *++argv;
  outputfilename  = *++argv;
  outputfilename2 = *++argv;
 
  inputfile3=efopen(inputfilename3,"r");
  fscanf(inputfile3,"%d",&numv);
  fscanf(inputfile3,"%d",&num);
  fclose(inputfile3);

  inputfile   =  efopen(inputfilename,"r");
  inputfile2  =  efopen(inputfilename2,"r");
  outputfile  = efopen(outputfilename,"w");
  outputfile2 = efopen(outputfilename2,"w");

  Lambda =(double *)ecalloc(numv*numv,sizeof(double));
  Sigma  =(double *)ecalloc(numv*numv,sizeof(double));
  Lambdaini=(double *)ecalloc(numv*numv,sizeof(double));
  Sigmaini =(double *)ecalloc(numv*numv,sizeof(double));
  KLdiv  =(double *)ecalloc(numv,sizeof(double));

  io_scanmatrix(inputfile,numv,numv,Sigma);
  io_scanmatrix(inputfile,numv,numv,Lambda);
  io_scanmatrix(inputfile2,numv,numv,Sigmaini);
  io_scanmatrix(inputfile2,numv,numv,Lambdaini);
  ssl_KLdiv(Lambda,Sigma,Lambdaini,Sigmaini,numv,KLdiv,topten,vtopten);
  outputKLdivtopten(outputfile,topten,vtopten,num);
  outputKLdiv(outputfile2,KLdiv,numv);

  fclose(inputfile);
  fclose(inputfile2);
  fclose(outputfile);
  fclose(outputfile2);
  free(Lambda);
  free(Lambdaini);
  free(Sigma);
  free(Sigmaini);

}

int outputKLdiv(FILE *outputfile,double *KLdiv,int numv){
  int i;

  for (i=0;i<numv;++i) {
    fprintf(outputfile,"%e ",KLdiv[i]);
  }
  fprintf(outputfile,"\n");    

}

int outputKLdivtopten(FILE *outputfile,int topten[10], double vtopten[10], int num){
  int i;

  for (i=0;i<num;++i)
    fprintf(outputfile,"%4d ",topten[i]);
  fprintf(outputfile,"\n     ");
  for (i=0;i<num;++i)
    fprintf(outputfile,"%e ",vtopten[i]);
  fprintf(outputfile,"\n");

}
