
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "SSL.h"
#include "EF.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k;
  int nums,numv,numite,numstotal,numslide;
  int topten[10],num;
  int inflag,outflag,outflag2;
  double *data,*datanorm;
  double *COVM,*Lambda,*Sigma,*Lambdap,*Sigmap,*Lambdaini,*Sigmaini;
  double rou,*KLdiv,vtopten[10],sumKLdiv;

  char *inputfilename,*inputfilename2;
  char *outputfilenamebase,*outputfilenamebase2,*outputfilename3,*outputfilename4,*outputfilename5,*outputfilename6;

  FILE *inputfile,*inputfile2;
  FILE *outputfile3,*outputfile4,*outputfile5,*outputfile6;
  
  char *line;
  size_t len=0;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;
  
  pi=acos(-1.0);
  
  struct option long_opt[] = {
    {"c",0,NULL,'c'},
    {"v",0,NULL,'v'},
    {"e",0,NULL,'e'},
    {"b",0,NULL,'b'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"h",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'c':
      break;
    case 'v':
      break;
    case 'e':
      break;
    case 'b':
      break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);  exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < ) {
    USAGE(progname);
    exit(1);
  }
  fscanf(inputfile2,"%d",&numstotal);
  fscanf(inputfile2,"%d",&nums);
  fscanf(inputfile2,"%d",&numslide);
  fscanf(inputfile2,"%d",&numv);
  fscanf(inputfile2,"%d",&num);
  fscanf(inputfile2,"%lf",&rou);
  fscanf(inputfile2,"%d",&outflag);
  fclose(inputfile2);

  = *argv;
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename \n",progname);
}


