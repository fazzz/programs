
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SSL.h"
#include "efunc.h"

int scnadata(FILE *inputfile,int numv,double *data);
int outputKLdivtopten(FILE *outputfile,int topten[10], double vtopten[10],double sumKLdiv);

int main(int argc, char *argv[]) {
  int i,j,k;
  int numv,topten[10];
  double *Lambda,*Lambda2,*KLdiv,*Sigma,*Sigmap,vtopten[10],sumKLdiv;

  char *inputfilename,*inputfilename2,*condfilename;
  char *outputfilename;
  FILE *inputfile,*inputfile2,*condfile;
  FILE *outputfile;

  if (argc < 3) {
    printf("USAGE: .ca_lssl.exe inputfilename(data) inputfilename2(data2) condfile(cond) outputfile \n");
    exit(1);
  }

  inputfilename   = *++argv;
  inputfilename2  = *++argv;
  condfilename    = *++argv;
  outputfilename  = *++argv;

  condfile=efopen(inputfilename2,"r");
  fscanf(inputfile2,"%d",&numv);
  fclose(inputfile2);

  inputfile   =  efopen(inputfilename,"r");
  inputfile2  =  efopen(inputfilename2,"r");
  outputfile  = efopen(outputfilename,"w");

  Lambda =(double *)ecalloc(numv*numv,sizeof(double));
  Lambda2=(double *)ecalloc(numv*numv,sizeof(double));
  KLdiv  =(double *)ecalloc(numv,sizeof(double));

  scnadata(inputfile,numv,datag);
  sumKLdiv=ssl_KLdiv(Lambda,Sigma,Lambdap,Sigmap,numv,KLdiv,topten,vtopten);
  output(outputfile,topten,vtopten,sumKLdiv);

  fclose(inputfile);
  fclose(outputfile3);
  free(Lambda);
  free(Lambda2);
  free(KLdiv);

  return 0;
}

int scnadata(FILE *inputfile,int numv,double *data){
  int i,j;

  for (i=0;i<numv;++i) {
    for (j=0;j<numv;++j) {
      fscanf(inputfile,"%lf",&f);
      data[i*numv+j]=f;
  }

  return 0;

}


int outputKLdivtopten(FILE *outputfile,int topten[10], double vtopten[10],double sumKLdiv){
  int i;

  for (i=0;i<10;++i) {
    fprintf(outputfile,"%4d ",topten[i]);
  }
  fprintf(outputfile,"\n     ");
  for (i=0;i<10;++i) {
    fprintf(outputfile,"%4.2lf ",vtopten[i]);
  }
  fprintf(outputfile,"\n");
  fprintf(outputfile,"%4.2lf \n",sumKLdiv);

}
