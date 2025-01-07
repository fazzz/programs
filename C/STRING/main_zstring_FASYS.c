
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "STRING.h"
#include "CSI.h"
#include "PT.h"
#include "FF.h"
#include "EF.h"
#include "IO.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int flag;

  int sl=0;

  int outinterval=1;

  int numpoint,numiteration;

  double dt=0.01;

  char *line;
  size_t len=0;

  char *progname;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *path,*path_evoluted;
  double *fe;

  double kd[2],n[2];

  struct potential e;
  struct force f;
  double *crd,*v;

  char *inifilename,*outputfilename,*outputfilename2;
  FILE *inifile,*outputfile,*outputfile2;

  kd[0]=2.5;
  kd[1]=2.5;
  n[0]=3.0;
  n[1]=3.0;

  progname=argv[0];

  while((c=getopt(argc,argv,"hs:t:o:k:l:n:m:"))!=-1) {
    switch(c) {
    case 's':
      sl=atoi(optarg);
      break;
    case 't':
      dt=atof(optarg);
      break;
    case 'o':
      outinterval=atoi(optarg);
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    case 'k':
      kd[0]=atof(optarg);
      break;
    case 'l':
      kd[1]=atof(optarg);
      break;
    case 'n':
      n[0]=atof(optarg);
      break;
    case 'm':
      n[1]=atof(optarg);
      break;
    default:
      USAGE(progname);
      exit(1);
    }
  }

  argc-=optind;
  argv+=optind;

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  numiteration=atoi(*argv);
  numpoint=atoi(*++argv);
  inifilename  = *++argv;
  outputfilename=*++argv;
  outputfilename2=*++argv;

  path=(double *)gcemalloc(sizeof(double)*5*3*numpoint);
  inifile=efopen(inifilename,"r");
  for (i=0;i<sl;++i) getline(&line,&len,inifile);
  io_scantrj2(inifile,5,numpoint,path);
  fclose(inifile);

  v=(double *)gcemalloc(sizeof(double)*numpoint);

  outputfile=efopen(outputfilename,"w"); 
  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<numiteration;++i) {
    z_string_FASYS(path,numpoint,dt,v,kd,n);
    if (i%outinterval==0) {
      io_outtrj2(outputfile,5,numpoint,path);
      fprintf(outputfile,"\n");
      
      for (j=0;j<numpoint;++j) fprintf(outputfile2,"%10.4lf\n",v[j]);
    }
  }
  fclose(outputfile);
  fclose(outputfile2);
}

int USAGE(char *progname) {
  printf("USAGE:%s \n",progname);
  printf("[-s sl] (dissmiss lines)    \n");
  printf("[-t dt]    \n");
  printf("[-o outinterval] \n");
  printf("[-k kd1] \n");
  printf("[-l kd2] \n");
  printf("[-m n1] \n");
  printf("[-n n2] \n");
  printf("[-h] help  \n");
  printf("%s  numiteration numpoint inifilename(path) outputfilename(path) outputfilename2(ene)\n",progname);
}

 
