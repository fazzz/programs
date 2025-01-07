#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "HIST.h"
#include "EF.h"

int scna1coldata(FILE *inputfile,int numf,int numi,int numcol,int col,double *data);

int main(int argc, char *argv[]) {
  int i,j;
  int numdata,numcol,col,normflag;
  int frame;
  double width;
  double max,min;
  double sum;
  double *data,*hist;

  char *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;

  if (argc < 7) {
    printf("USAGE: ./%s normflag[o(n)(of)fa(ve)(rec)t] inputfilename(data) outputfilename(hist) width numdata numcol col \n",argv[0]);
    exit(1);
  }
  normflag=(*++argv)[0];
  if (normflag != 'o' && normflag != 'f' && normflag != 'a' && normflag != 't') {
    printf("flag error: must be o(n)  or (of)f or a(ve) or (rec)t "); exit(1);
  }
  inputfilename   = *++argv;
  outputfilename  = *++argv;
  width  = atof(*++argv);
  numdata= atoi(*++argv);
  numcol = atoi(*++argv);
  if (argc < 7)
    col=1;
  else
    col = atoi(*++argv);

  inputfile   =  efopen(inputfilename,"r");
  outputfile  = efopen(outputfilename,"w");

  data=(double *)gcemalloc(sizeof(double)*numdata);
  scna1coldata(inputfile,numdata,1,numcol,col,data);
  hist=hist_mkhist(data,width,numdata,&max,&min,&frame,normflag);

  sum=0.0;
  for (i=0;i<=frame;++i) sum+=hist[i];

  if (normflag=='f')
    for (i=0;i<=frame;++i) fprintf(outputfile,"%e %e \n",width*(i+0.5)+min,hist[i]);
  else if (normflag=='o')
    for (i=0;i<=frame;++i) fprintf(outputfile,"%e %e \n",width*(i+0.5)+min,hist[i]/sum/width);
  else if (normflag=='a')
    for (i=0;i<=frame;++i) fprintf(outputfile,"%e %e \n",width*(i+0.5)+min,hist[i]/width);
  else {
    for (i=0;i<=frame;++i) {
      fprintf(outputfile,"%e %e \n",width*i+min,0.0);
      fprintf(outputfile,"%e %e \n",width*(i+1.0)+min,hist[i]/width);
    }
  }
  fclose(inputfile);
  fclose(outputfile);
}

int scna1coldata(FILE *inputfile,int numf,int numi,int numcol,int col,double *data){
  int i,j;
  double temp;
  char *line;
  size_t len=0;

  for (i=0;i<numi;++i)
    getline(&line,&len,inputfile);  
  for (i=0;i<numf-numi+1;++i) {
    for (j=0;j<numcol;++j) {
      if (j==col-1)
	fscanf(inputfile,"%lf",&data[i]);
      else 
	fscanf(inputfile,"%lf",&temp);
    }
  }

  return 0;
}

