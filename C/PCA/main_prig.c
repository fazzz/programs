#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PT.h"
#include "PCA.h"
#include "IO.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j;
  int flag=0;
  int numstep,numterm;
  double *intene,*cov,*eigenval,*pepca,tp,sum,sum2;

  char *inputfilename1,*inputfilename2,*inputfilename3,*inputfilename4,*outputfilename;
  FILE *inputfile1,*inputfile2,*inputfile3,*inputfile4,*outputfile;
  
  if (argc < 5) {
    printf("USAGE: %s flag(n or d or b) inputfilename1(data) inputfilename2(cov) inputfilename3(pt) inputfilename4(cond.txt) outputfilename \n",argv[0]);
    exit(1);
  }
  flag=(*++argv)[0];
  if (flag != 'n' && flag != 'd' && flag != 'b') {
    printf("flag error: must be n  or d or b ");
    exit(1);
  }
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  inputfilename3 = *++argv;
  inputfilename4 = *++argv;
  outputfilename = *++argv;
  
  inputfile4=efopen(inputfilename4,"r");
  fscanf(inputfile4,"%d",&numterm);
  fscanf(inputfile4,"%lf",&tp);
  fscanf(inputfile4,"%d",&numstep);
  fclose(inputfile4);
  inputfile4=efopen(inputfilename4,"r");
  readParmtop(inputfile4);
  fclose(inputfile4);

  intene=(double *)gcemalloc(sizeof(double)*numstep*numterm);
  pepca=(double *)gcemalloc(sizeof(double)*numstep*numterm);
  cov=(double *)gcemalloc(sizeof(double)*numterm*numterm);
  for (i=0;i<numstep;++i) {
    for (j=0;j<numterm;++j){
      intene[i*numterm+j]=0.0;
    }  
  }
  inputfile1=efopen(inputfilename1,"r");
  io_scantimeseries(inputfile1,numstep,numterm,intene,'c');
  fclose(inputfile1);

  for (i=0;i<numterm;++i)
    for (j=0;j<numterm;++j)
      cov[i*numterm+j]=0.0;
  inputfile2=efopen(inputfilename2,"r");
  for (i=0;i<numterm;++i)
    for (j=0;j<numterm;++j)
      fscanf(inputfile2,"%lf",&cov[i*numterm+j]);
  fclose(inputfile2);

  pepca_norm(intene,numstep,numterm,tp);
  pepca_proj(intene,pepca,cov,numstep,numterm);
  
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    fprintf(outputfile,"%d ",i);
    for (j=0;j<numterm;++j) {
      fprintf(outputfile,"%lf ",pepca[i*numterm+j]);
    }
    fprintf(outputfile,"\n");
  }
  fclose(outputfile);

  return 0;
}
