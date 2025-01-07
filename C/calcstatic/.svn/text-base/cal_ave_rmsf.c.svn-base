
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "EF.h"
#include "mymath.h"

int io_scanodtimeseries(FILE *inputfile,int numcolum,int numspecificcolum,int *numspecificcolumlist, int numintialraw, int numfinalrow, double *data);

int main(int argc, char *argv[]) {
  int i,j;
  int numcol,numv,numi,numf,numstep;
  int *numvlist;
  double *ave,*rmsf;
  double *data,*temp;
  
  char *inputfilename1,*inputfilename2,*outputfilename;
  FILE *inputfile1,*inputfile2,*outputfile;
  
  if (argc < 3) {
    printf("USAGE: ave_rmsf inputfilename1(data) inputfilename2(cond)\n");
    printf("cond: numcol numv numvlist(......) numi numf\n");
    exit(1);
  }
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  if (argc==4)
    outputfilename = *++argv;
  
  inputfile2=efopen(inputfilename2,"r");
  fscanf(inputfile2,"%d",&numcol);
  fscanf(inputfile2,"%d",&numv);
  numvlist=(int *)ecalloc(sizeof(int),numv);
  for (i=0;i<numv;++i) {
    fscanf(inputfile2,"%d",&numvlist[i]);--numvlist[i];
  }
  fscanf(inputfile2,"%d",&numi);--numi;
  fscanf(inputfile2,"%d",&numf);--numf;
  fclose(inputfile2);
  numstep=numf-numi+1;
  data=(double *)ecalloc(sizeof(double),numstep*numv);
  
  inputfile1=efopen(inputfilename1,"r");
  io_scanodtimeseries(inputfile1,numcol,numv,numvlist,numi,numf,data);
  fclose(inputfile1);

  ave =(double *)ecalloc(sizeof(double),numv);
  rmsf=(double *)ecalloc(sizeof(double),numv);
  for (i=0;i<numv;++i) {
    temp=(double *)ecalloc(sizeof(double),numstep);
    for (j=0;j<numstep;++j)
      temp[j]=data[i*numstep+j];
    ave[i] =calc_ave(numstep,temp);
    rmsf[i]=calc_flu(numstep,temp);
    free(temp);
  }

  if (argc>3) {
    outputfile=efopen(outputfilename,"w");
    fprintf(outputfile,"%s ",inputfilename1);
    for (i=0;i<numv;++i) {
      fprintf(outputfile,"%d ",numvlist[i]+1);
    }
    fprintf(outputfile,"\nave ");
    for (i=0;i<numv;++i) {
      fprintf(outputfile,"%lf ",ave[i]);
     }
    fprintf(outputfile,"\nrmsf ");
    for (i=0;i<numv;++i) {
      fprintf(outputfile,"%lf ",rmsf[i]);
    }
    fprintf(outputfile,"\n ");
    fclose(outputfile);
  }
  else {
    printf("%s ",inputfilename1);
    for (i=0;i<numv;++i) {
      printf("%d ",numvlist[i]+1);
    }
    printf("\nave ");
    for (i=0;i<numv;++i) {
      printf("%lf ",ave[i]);
     }
    printf("\nrmsf ");
    for (i=0;i<numv;++i) {
      printf("%lf ",rmsf[i]);
    }
    printf("\n ");
  }
  free(data);
  free(numvlist);
  free(ave);
  free(rmsf);
  
  return 0;
}



int io_scanodtimeseries(FILE *inputfile,int numcol,int numv,int *numvlist, int numi, int numf, double *data){
  int i,j,k,l=0;
  double temp;
  char *line;
  size_t len=0;

  for (i=0;i<numi;++i)
    getline(&line,&len,inputfile);  
  for (i=0;i<numf-numi+1;++i) {
    for (j=0;j<numcol;++j) {
      for (k=0;k<numv;++k) {
	if (j==numvlist[k]) {
	  fscanf(inputfile,"%lf",&data[k*(numf-numi+1)+i]);
	  ++l;
	  break;
	}
      }
      if (k==numv)
	fscanf(inputfile,"%lf",&temp);
    }
  }
  
  return 0;
}
