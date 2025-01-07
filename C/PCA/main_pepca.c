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
  double *intene,*cov,*eigenval,*pepca,tp,sum,sum2,*ave,*var;

  char *inputfilename1,*inputfilename2,*inputfilename3,*outputfilename,*outputfilename2,*outputfilename3;
  FILE *inputfile1,*inputfile2,*inputfile3,*outputfile,*outputfile2,*outputfile3,*outputfile4;
  
  if (argc < 6) {
    printf("USAGE: %s flag(n or d or b) inputfilename1(data) inputfilename2(pt) inputfilename3(cond.txt) outputfilename outputfilename2 outputfilename3\n",argv[0]);
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
  outputfilename = *++argv;
  outputfilename2= *++argv;
  outputfilename3= *++argv;
  
  inputfile3=efopen(inputfilename3,"r");
  fscanf(inputfile3,"%d",&numterm);
  fscanf(inputfile3,"%lf",&tp);
  fscanf(inputfile3,"%d",&numstep);
  fclose(inputfile3);
  inputfile2=efopen(inputfilename2,"r");
  readParmtop(inputfile2);
  fclose(inputfile2);
  //  numterm=(AP.NATOM*AP.NATOM)-AP.NEXT;
  intene=(double *)gcemalloc(sizeof(double)*numstep*numterm);
  pepca=(double *)gcemalloc(sizeof(double)*numstep*numterm);
  cov=(double *)gcemalloc(sizeof(double)*numterm*numterm);
  eigenval=(double *)gcemalloc(sizeof(double)*numterm);
  ave=(double *)gcemalloc(sizeof(double)*numterm);
  var=(double *)gcemalloc(sizeof(double)*numterm);
  for (i=0;i<numstep;++i) {
    for (j=0;j<numterm;++j){
      intene[i*numterm+j]=0.0;
    }  
  }
  for (i=0;i<numterm;++i)
    for (j=0;j<numterm;++j)
      cov[i*numterm+j]=0.0;
  inputfile1=efopen(inputfilename1,"r");
  io_scantimeseries(inputfile1,numstep,numterm,intene,'c');
  fclose(inputfile1);

  pepca_norm(intene,numstep,numterm,tp);
  pepca_covm(intene,numstep,numterm,cov);
  pepca_diag(cov,eigenval,numterm);
  pepca_proj(intene,pepca,cov,numstep,numterm);
  pepca_avevar(pepca,numstep,numterm,ave,var);
  
  sum=0.0;sum2=0;
  for (i=0;i<numterm;++i)
    sum+=eigenval[i];
  outputfile=efopen(outputfilename,"w");
  fprintf(outputfile,"n eigenvalue(sigma^2) %% \n");
  for (i=0;i<numterm;++i) {
    sum2+=eigenval[i];
    fprintf(outputfile,"%d %12.8lf %12.8lf \n",i,eigenval[i],sum2/sum*100.0);
  }
  fclose(outputfile);

  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<numterm;++i) {
    for (j=0;j<numterm;++j)
      fprintf(outputfile2,"%lf ",cov[i*numterm+j]);
    fprintf(outputfile2,"\n");
  }
  fclose(outputfile2);

  outputfile3=efopen(outputfilename3,"w");
  for (i=0;i<numstep;++i) {
    fprintf(outputfile3,"%d ",i);
    for (j=0;j<numterm;++j) {
      fprintf(outputfile3,"%lf ",pepca[i*numterm+j]);
    }
    fprintf(outputfile3,"\n");
  }
  fclose(outputfile3);

  outputfile4=efopen("log_av_pepca.txt","w");
  for (i=0;i<numterm;++i) {
    fprintf(outputfile4,"%lf ",ave[i]);
    fprintf(outputfile4,"%lf \n",var[i]);
  }
  fclose(outputfile4);


  return 0;
}
