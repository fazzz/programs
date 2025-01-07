
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EF.h"
#include "RAND.h"

int USAGE(char *progname);

double Box_Muller2(int i,double mean,double vari);

int main(int argc, char *argv[]) {
  int i,j,k;

  int N;

  double f;

  int K;
  double *nyu_k,*Sigma_k, *pi_k;

  double pi;

  char *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;
  
  char *line;
  size_t len=0;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;
  
  pi=acos(-1.0);
  
  struct option long_opt[] = {
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"h",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);  exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  K = atoi(*argv);
  N = atoi(*++argv);
  inputfilename  = *++argv;
  outputfilename = *++argv;

  nyu_k=(double *)gcemalloc(sizeof(double)*K);
  Sigma_k=(double *)gcemalloc(sizeof(double)*K);

  inputfile=efopen(inputfilename,"r");
  for (i=0;i<K;++i) {
    fscanf(inputfile,"%lf",&nyu_k[i]);
    fscanf(inputfile,"%lf",&Sigma_k[i]);
  }
  fclose(inputfile);

  N=(int)(N/K);
  
  outputfile=efopen(outputfilename,"w");
  k=0;
  for (i=0;i<N;++i) {
    for (j=0;j<K;++j) {
      ++k;
      f=Box_Muller2(i,nyu_k[j],Sigma_k[j]);
      fprintf(outputfile,"%3d %10.8lf \n",k,f);
    }
  }
  fclose(outputfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] N K inputfilename outputfilename \n",progname);
}


double Box_Muller2(int i,double mean,double vari) {
  double u1,u2,pi;
  double c1,c2;
  //  FILE *logfile;

  u1=genrand_real2();
  u2=genrand_real2();
  pi=acos(-1);

  /**************************************************************************/
  /* logfile=efopen("log.txt","a");					    */
  /* c1=sqrt(2.0)*cos(2.0*pi*u2);					    */
  /* c1=log(2.0);							    */
  /* c1=sqrt(-2.0*log(u1))*cos(2.0*pi*u2);				    */
  /* c2=sqrt(-2.0*log(u1))*sin(2.0*pi*u2);				    */
  /* fprintf(logfile,"%d %12.8lf %12.8lf %12.8lf %12.8lf\n",i,u1,u2,c1,c2); */
  /* fclose(logfile);							    */
  /**************************************************************************/


  if (i%2) {
    c1=mean+sqrt(-2.0*log(u1))*cos(2.0*pi*u2)*sqrt(vari);
    return c1;
  }
  else {
    c2=mean+sqrt(-2.0*log(u1))*sin(2.0*pi*u2)*sqrt(vari);
    return c2;
  }
}
