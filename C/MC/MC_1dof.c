
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EF.h"

#include "MC.h"
#include "RAND.h"

#define ON 1
#define OFF 0

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j;
  int numstep=100000,interval=1;
  
  int rstflag=OFF;
  double beta=1.0;
  double x=0.0,x_trial,dx=10.0;
  double V,V_trial,delta;

  char *outputfilename,*trjfilename;
  char *rstfilename,*outrstfilename;
  FILE *outputfile,*trjfile;
  FILE *rstfile,*outrstfile;

  char *progname;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  int opt_idx=1;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  struct option long_opt[] = {
    {"dx",1,NULL,'d'},
    {"int",1,NULL,'i'},
    {"nums",1,NULL,'s'},
    {"rst",1,NULL,'r'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hs:d:i:t:r:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 't':
      beta=atof(optarg);
      break;
    case 's':
      numstep=atoi(optarg);
      break;
    case 'd':
      dx=atof(optarg);
      break;
    case 'i':
      interval=atoi(optarg);
      break;
    case 'r':
      rstflag=ON;
      rstfilename=optarg;
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  outputfilename    = *argv;
  trjfilename       = *++argv;
  outrstfilename    = *++argv;

  outputfile=efopen(outputfilename,"w");
  trjfile=efopen(trjfilename,"w");

  if (rstflag==ON) {
    rstfile=efopen(rstfilename,"r");
    fscanf(rstfile,"%lf",&x);
    fclose(rstfile);
  }

  V=0.5*x*x;

  for (i=0;i<numstep;++i) {
    x_trial=x+dx*(genrand_real2()-0.5);

    V_trial=0.5*(x_trial*x_trial);

    delta=V_trial-V;
    
    if((c=Metropolis(beta*delta))==1) {
      x=x_trial;
      V=V_trial;
    }

    if (i%interval==0) {
      fprintf(outputfile,"%10.8lf %d\n",V,c);
      fprintf(trjfile,"%10.8lf %10.8lf \n",x);
    }
  }
  fclose(outputfile);
  fclose(trjfile);

  outrstfile=efopen(outrstfilename,"w");
  fprintf(outrstfile,"%lf %lf\n",x);
  fclose(outrstfile);


  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] outputfilename trjfilename outrstfilename\n",progname);
}


