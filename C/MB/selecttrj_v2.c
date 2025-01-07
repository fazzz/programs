
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "IO.h"
#include "PT.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j;
  int numstep,numatom,numdihed;
  int num;
  int totalnum;
  int **atomdihedpairs;
  double max,min,interval;
  double pca,*data;
  char *line;
  size_t len=0;
  
  char *inputfilename1,*inputfilename2,*inputfilename3,*inputfilename4,*outputfilenamebase,outputfilename[100];
  FILE *inputfile1,*inputfile2,*inputfile3,*inputfile4,*outputfile;
  
  if (argc < 6) {
    printf("USAGE: %s interval(kcal/mol) inputfilename1(dihed) inputfilename2(pepca) cond inputfilename4(parm) outputfilenamebase \n",argv[0]);
    exit(1);
  }
  /**********************/
  /* max=atof(*++argv); */
  /* min=atof(*++argv); */
  /**********************/
  interval=atof(*++argv);
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  inputfilename3 = *++argv;
  inputfilename4 = *++argv;
  outputfilenamebase = *++argv;

  inputfile3=efopen(inputfilename3,"r");
  fscanf(inputfile3,"%d",&numstep);
  fscanf(inputfile3,"%d",&numdihed);
  fclose(inputfile3);

  inputfile2=efopen(inputfilename2,"r");
  fscanf(inputfile2,"%lf",&pca);
  min=pca;max=pca;
  for (i=1;i<numstep;++i) {
    fscanf(inputfile2,"%lf",&pca);
    if (pca < min) min=pca;
    if (pca > max) max=pca;
  }
  fclose(inputfile2);

  printf("max=%lf\n",max);
  printf("min=%lf\n",min);
  
  //  interval = (max-min)/totalnum;

  inputfile4=efopen(inputfilename4,"r");
  readParmtop(inputfile4);
  fclose(inputfile4);

  data=(double *)gcemalloc(sizeof(double)*(numdihed+1));
  
  inputfile1=efopen(inputfilename1,"r");
  getline(&line,&len,inputfile1);
  inputfile2=efopen(inputfilename2,"r");

  for (i=0;i<=numstep;++i) {
    io_scan_data(inputfile1,numdihed+1,data,'x');
    fscanf(inputfile2,"%lf",&pca);

    num=(int)((pca-min)/interval);
    printf("%d\n",num);
    sprintf(outputfilename,"%s_%d",outputfilenamebase,num);
    outputfile=efopen(outputfilename,"a");
    io_outputdata(outputfile,numdihed+1,data);
    fclose(outputfile);
  }
  
  fclose(inputfile1);
  fclose(inputfile2);
  
  return 0;
}

