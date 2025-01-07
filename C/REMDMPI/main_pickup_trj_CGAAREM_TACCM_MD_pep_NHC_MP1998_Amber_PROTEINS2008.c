#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EF.h"

#include "REMDCGAA_TACCM_calc_uene_Amber_PROTEINS2008.h"

#include "REMDCGAA_TACCM_MPI_2_Amber_PROTEINS2008.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,m,n,d;
  int numd=1;
  int index;
  int *numbering;

  int numstep=1000,interval=100;
  int numEX=1,numRE=2;

  int equflag=OFF;
  int numstepequ=0,intervalequ=1;

  double *trj;

  int **sereis;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double pi;

  char **trjfilename;
  char *inputfilename;
  FILE *inputfile;

  char *parametertrjfilename;
  FILE *parametertrjfile;

  char *outputfilenamebase,outputfilename[2000];

  FILE **trjfile,**outputfile;

  char *progname;

  int opt_idx=1;

  pi=acos(-1.0);

  struct option long_opt[] = {
    {"numRE",1,NULL,'N'}, {"numEX",1,NULL,'e'},
    {"nums",1,NULL,'S'}, {"int",1,NULL,'i'},
    {"equ",1,NULL,'E'}, {"intequ",1,NULL,'I'},
    {"numd",1,NULL,'n'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hnN:e:S:i:n:E:I:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'e':
      numEX=atoi(optarg); break;
    case 'N':
      numRE=atoi(optarg); break;
    case 'S':
      numstep=atoi(optarg); break;
    case 'i':
      interval=atoi(optarg); break;
    case 'E':
      equflag=ON;
      numstepequ=atoi(optarg);  break;
    case 'I':
      intervalequ=atoi(optarg);  break;
    case 'n':
      numd=atoi(optarg);  break;
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
  inputfilename        = *argv;
  parametertrjfilename = *++argv;
  outputfilenamebase   = *++argv;

  trj=(double *)gcemalloc(sizeof(double)*numd);

  trjfilename=(char **)gcemalloc(sizeof(char *)*numRE);
  for (i=0;i<numRE;++i) trjfilename[i]=(char *)gcemalloc(sizeof(char)*1000);
  trjfile=(FILE **)gcemalloc(sizeof(trjfile)*numRE);

  inputfile=efopen(inputfilename,"r");
  CGAAREMDreadInputs_pickup_trj(inputfile,numRE,trjfilename,trjfile);
  fclose(inputfile);

  sereis=(int **)gcemalloc(sizeof(int *)*numRE);
  for (i=0;i<numRE;++i) sereis[i]=(int *)gcemalloc(sizeof(int)*numEX);
  numbering=(int *)gcemalloc(sizeof(int)*numRE);
  for (i=0;i<numRE;++i) numbering[i]=0;

  parametertrjfile=efopen(parametertrjfilename,"r");
  for (i=0;i<numRE;++i)
    for (j=0;j<numEX;++j) 
      fscanf(parametertrjfile,"%d",&(sereis[i][j]));
  fclose(parametertrjfile);

  outputfile=(FILE **)gcemalloc(sizeof(FILE *)*numRE);
  for (i=0;i<numRE;++i) {
    sprintf(outputfilename,"%s_%d\0",outputfilenamebase,i+1);
    outputfile[i]=efopen(outputfilename,"w");
  }

  numstep=numstep/interval;
  numstepequ=numstepequ/intervalequ;
  for (i=0;i<numRE;++i) {
    if (equflag==ON) 
      for (j=0;j<numstepequ;++j) 
	for (k=0;k<numd;++k) fscanf(trjfile[i],"%lf",&trj[k]);

    for (j=0;j<numEX;++j) {
      index=sereis[i][j];

      for (k=0;k<numstep;++k) {
	for (m=0;m<numd;++m) fscanf(trjfile[i],"%lf",&trj[m]);

	numbering[index]+=1;
	for (m=0;m<numd;++m)  fprintf(outputfile[index],"%e ",trj[m]);
	fprintf(outputfile[index],"\n");
      }
    }
  }
  for (i=0;i<numRE;++i) fclose(outputfile[i]);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parametertrjfilename parmfilename TACCMfilename outputfilename\n",progname);
}

