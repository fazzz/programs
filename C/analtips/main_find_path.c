#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EF.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,n,m,n_o,m_o;

  int numbinx,numbiny,numpath;

  double norm=1.0;

  double num,xmin,ymin,xmax,ymax,widthx,widthy;

  double **x,**y,**pmf,**peo,**path;
  
  char *inputfilename,*pathfilename,*outputfilename;
  FILE *inputfile,*pathfile,*outputfile;

  char *line;
  size_t len=0;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;
  
  struct option long_opt[] = {
    {"n",1,NULL,'n'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"hn:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'n':
      norm=atof(optarg);
      break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);   exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }
  numbinx = atoi(*argv);
  numbiny = atoi(*++argv);
  numpath = atoi(*++argv);
  pathfilename  = *++argv;
  inputfilename = *++argv;
  outputfilename = *++argv;

  x=(double **)gcemalloc(sizeof(double *)*numbinx);
  for (i=0;i<numbinx;++i) x[i]=(double *)gcemalloc(sizeof(double)*numbiny); 
  y=(double **)gcemalloc(sizeof(double *)*numbinx);
  for (i=0;i<numbinx;++i) y[i]=(double *)gcemalloc(sizeof(double)*numbiny); 
  pmf=(double **)gcemalloc(sizeof(double *)*numbinx);
  for (i=0;i<numbinx;++i) pmf[i]=(double *)gcemalloc(sizeof(double)*numbiny); 
  peo=(double **)gcemalloc(sizeof(double *)*numbinx);
  for (i=0;i<numbinx;++i) peo[i]=(double *)gcemalloc(sizeof(double)*numbiny); 

  inputfile=efopen(inputfilename,"r");
  for (i=0;i<numbinx;++i) {
    for (j=0;j<numbiny;++j) {
      fscanf(inputfile,"%lf",&x[i][j]);
      fscanf(inputfile,"%lf",&y[i][j]);
      fscanf(inputfile,"%lf",&pmf[i][j]);
      fscanf(inputfile,"%lf",&peo[i][j]);
    }
  }
  fclose(inputfile);

  xmin=x[0][0];
  ymin=y[0][0];
  xmax=x[numbinx-1][numbiny-1];
  ymax=y[numbinx-1][numbiny-1];

  widthx=(xmax-xmin)/numbinx;
  widthy=(ymax-ymin)/numbiny;

  path=(double **)gcemalloc(sizeof(double *)*numpath);
  for (i=0;i<numpath;++i) path[i]=(double *)gcemalloc(sizeof(double)*2); 

  pathfile=efopen(pathfilename,"r");
  for (i=0;i<numpath;++i) {
    for (j=0;j<2;++j) {
      fscanf(pathfile,"%lf",&path[i][j]);
    }
  }
  fclose(pathfile);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numpath;++i) {
    n=(int)((path[i][0]-xmin)/widthx);
    m=(int)((path[i][1]-ymin)/widthy);

    num=((double)(i))/numpath*norm;

    if (n!=n_o || m!=m_o) {
      fprintf(outputfile,"%12.4lf %12.4lf %12.4lf %12.4lf %12.4lf\n",num,x[n][m],y[n][m],pmf[n][m],peo[n][m]);
    }

    n_o=n;
    m_o=m;
  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] numbinx numbiny numpath pathfilename inputfilename outputfilename\n",progname);
}


