
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PMF.h"
#include "PT.h"
#include "EF.h"

double *io_scandcoldata2(FILE *inputfile,int numi,int numcol,int xcol,int ycol,int *numstep,double *data);

int main(int argc, char *argv[]) {
  int i,j;

  int numstep;
  double *data;
  int numi,numcol,xcol,ycol;
  double maxx,maxy,minx,miny;
  int framex,framey;
  double width;
  double *pmf;
  
  char *inputfilename,*outputfilename,*c;
  FILE *inputfile1,*outputfile;
  
  if (argc < 8) {
    printf("USAGE: %s  width numi numcol xcol ycol inputfilenames(2d_data) outputfilename(pmf2d) \n",argv[0]);
    exit(1);
  }
  width=atof(*++argv);
  numi=atoi(*++argv);
  numcol=atoi(*++argv);
  xcol=atoi(*++argv);
  ycol=atoi(*++argv);
  inputfilename  = *++argv;
  outputfilename = *++argv;
  numstep=0;
  data=(double *)gcemalloc(sizeof(double)*2);
  while ((c=*++argv) != NULL) {
    inputfilename=c;
    inputfile1=efopen(inputfilename,"r");
    data=io_scandcoldata2(inputfile1,numi,numcol,xcol,ycol,&numstep,data);
    fclose(inputfile1);
  }  

  pmf=pmf_2dmap(data,numstep,width,&maxx,&maxy,&minx,&miny,&framex,&framey);
  
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<framex;++i) {
    for (j=0;j<framey;++j)
      fprintf(outputfile,"%e %e %e\n",width*i+minx,width*j+miny,pmf[i*framey+j]);
    fprintf(outputfile,"\n");
  }
  fclose(outputfile);
  
  return 0;
}

double *io_scandcoldata2(FILE *inputfile,int numi,int numcol,int xcol,int ycol,int *numstep,double *data){
  int i,j,k;
  double f;

  char *line;
  size_t len=0;

  for (i=0;i<numi;++i)
    getline(&line,&len,inputfile);

  for (i=(*numstep);;++i) {
    for (j=0;j<numcol;++j) {
      if (fscanf(inputfile,"%lf",&f)!=-1) {
	if (j==xcol-1) {
	  data=(double *)gcerealloc(data,sizeof(double)*(i+1)*2);
	  data[i*2]=f;
	}
	else if (j==ycol-1) {
	  data[i*2+1]=f;
	}
      }
      else {
	*numstep=i-1;
	return data;
      }
    }
  }
}
