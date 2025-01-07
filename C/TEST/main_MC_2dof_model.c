
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
#include "BOXMULL.h"

#define ON 1
#define OFF 0

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j;
  int numstep=100000,interval=1000;
  
  int rstflag=OFF;
  double x=0.0,y=0.0,x_trial,y_trial,dx=0.001;
  double V,V_trial,delta;

  double beta=1,kb=1.98723e-3;

  double D=5.0,a=1.0,K=1.0,lamda=2.878;

  char *rstfilename;
  char *outputfilename,*trjfilename,*outrstfilename;
  char x_trjfilename[100],y_trjfilename[100];
  FILE *outputfile,*outrstfile;
  FILE *x_trjfile,*y_trjfile;
  FILE *rstfile;

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
    {"temp",1,NULL,'t'},
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

  beta=1.0/(kb*beta);

  sprintf(x_trjfilename,"%s_x",trjfilename);
  sprintf(y_trjfilename,"%s_y",trjfilename);

  if (rstflag==ON) {
    rstfile=efopen(rstfilename,"r");
    fscanf(rstfile,"%lf %lf",&x,&y);
    fclose(rstfile);
  }

  outputfile=efopen(outputfilename,"w");  
  x_trjfile=efopen(x_trjfilename,"w");
  y_trjfile=efopen(y_trjfilename,"w");
  fprintf(outputfile,"%10.8lf %d\n",V,c);
  fprintf(x_trjfile,"%10.8lf \n",x);
  fprintf(y_trjfile,"%10.8lf \n",y);

  V=D*(x*x-a*a)*(x*x-a*a)+0.5*K*y*y+lamda*x*y;
  //  V=0.5*D*(x-a)*(x-a)+0.5*K*y*y;

  for (i=0;i<numstep;++i) {
    x_trial=x+dx*(genrand_real2()-0.5);
    y_trial=y+dx*(genrand_real2()-0.5);

    V_trial=D*(x_trial*x_trial-a*a)*(x_trial*x_trial-a*a)
      +0.5*K*y_trial*y_trial+lamda*x_trial*y_trial;
    //    V_trial=0.5*D*(x_trial-a)*(x_trial-a)+0.5*K*y_trial*y_trial;

    delta=V_trial-V;
    
    if((c=Metropolis(beta*delta))==1) {
      x=x_trial;
      y=y_trial;
      V=V_trial;
    }

    if (i%interval==0) {
      fprintf(outputfile,"%10.8lf %d\n",V,c);
      fprintf(x_trjfile,"%10.8lf \n",x);
      fprintf(y_trjfile,"%10.8lf \n",y);
    }
  }
  fclose(outputfile);
  fclose(x_trjfile);
  fclose(y_trjfile);

  outrstfile=efopen(outrstfilename,"w");
  fprintf(outrstfile,"%10.8lf %10.8lf\n",x,y);
  fclose(outrstfile);


  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] outputfilename trjfilename outrstfilename\n",progname);
}


