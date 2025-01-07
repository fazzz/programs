
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#define ON 0
#define OFF 1

#define INTEGER 0
#define DECIMAL 1
#define SPACE   2
#define WRITE   3

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,n,m;

  int flag=OFF;

  int status,nf,nl,numstep;

  double d,f,v,sign;

  double x1,x2,y1,y2;

  double *datain,*dataout;
  double cos,sin,a;

  double pi;
  
  char *line;
  size_t len=0;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;
  
  char *progname;
  int opt_idx=1;
  
  pi=acos(-1.0);
  
  struct option long_opt[] = {
    {"f",0,NULL,'f'},
    {"x1",1,NULL,'1'},
    {"y1",1,NULL,'2'},
    {"x2",1,NULL,'3'},
    {"y2",1,NULL,'4'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"hf1:2:3:4:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'f':
      flag=ON;  
      break;
    case '1':
      x1=atof(optarg);  
      break;
    case '2':
      y1=atof(optarg);  
      break;
    case '3':
      x2=atof(optarg);  
      break;
    case '4':
      y2=atof(optarg);  
      break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);   exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  numstep = atoi(*argv);
  inputfilename=  *++argv;
  outputfilename= *++argv;

  inputfile=efopen(inputfilename,"r");

  datain=(double *)gcemalloc(sizeof(double)*numstep*2);

  if (flag==ON)
    getline(&line,&len,inputfile);
  for (i=0;i<numstep;++i) {
    if (flag==ON)
      fscanf(inputfile,"%d ",&m);
    fscanf(inputfile,"%lf %lf",&datain[i*2],&datain[i*2+1]);
  }

  fclose(inputfile);

  a=(y2-y1)/(x2-x1);

  cos=a/(sqrt(1.0+a*a));
  sin=sqrt(1.0-cos*cos);

  dataout=(double *)gcemalloc(sizeof(double)*numstep*2);

  for (i=0;i<numstep;++i) {
    dataout[i*2]=datain[i*2]*cos-datain[i*2+1]*sin;
    dataout[i*2+1]=datain[i*2]*sin+datain[i*2+1]*cos;
  }  

  outputfile=efopen(outputfilename,"w");
  if (flag==ON)
    fprintf(outputfile,"numstep 1 2\n ");
  for (i=0;i<numstep;++i) {
    if (flag==ON)
      fprintf(outputfile,"%d ",i);
    fprintf(outputfile,"%7.6e %7.6e\n",dataout[i*2],dataout[i*2+1]);
  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename \n",progname);
}
