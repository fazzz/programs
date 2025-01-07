
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>

#include "fftw3.h"

#include "IO.h"
#include "PT.h"
#include "EF.h"
#include "SPE.h"
#include "const.h"

#include "netcdf_mine.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,n;
  int flag='c';

  int numdihed=0;
  int num_of_split=0,num_step_ini=1,num_step_fin=0,numstep;

  double omg;
  double cv=2.999792e-2;
  double kb=1.98723e-3*4.18407*100.0;

  double deltat=0.001,pi,kbT,temp=300;
  double *spe,*trj,*speave,*sumspe;

  char *inputfilename1,*outputfilename,*progname;
  FILE *inputfile1,*outputfile;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *line;
  size_t len=0;

  while((c=getopt(argc,argv,"hcvn:d:s:i:f:t:"))!=-1) {
    switch(c) {
    case 'v':
      flag='v';
      break;
    case 'c':
      flag='c';
      break;
    case 'n':
      numdihed=atoi(optarg);
      break;
    case 'd':
      deltat=atof(optarg);
      break;
    case 's':
      num_of_split=atoi(optarg);
      break;
    case 'i':
      num_step_ini=atoi(optarg);
      break;
    case 'f':
      num_step_fin=atoi(optarg);
      break;
    case 't':
      temp=atof(optarg);
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  pi=acos(-1.0);

  if (numdihed==0) {
    printf("numdihed must not be 0\n");
    exit(1);
  }
  else if (num_of_split==0) {
    printf("num must not be 0\n");
    exit(1);
  } 
  else if (num_step_fin==0) {
    printf("num must not be 0\n");
    exit(1);
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;
  
  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  inputfilename1 = *argv;
  outputfilename = *++argv;
  
  kbT=kb*temp;

  numstep=(num_step_fin-num_step_ini)/num_of_split;
  inputfile1=efopen(inputfilename1,"r");
  io_dismissdata(inputfile1,num_step_ini-1,numdihed);
  speave=(double *)gcemalloc(sizeof(double)*numstep*numdihed);
  for (i=0;i<num_of_split;++i) {
    spe=(double *)gcemalloc(sizeof(double)*numstep*numdihed);
    trj=(double *)gcemalloc(sizeof(double)*numstep*numdihed);
    io_scandtraj(inputfile1,numstep,numdihed,trj);
    CSPd_decom(numdihed,numstep,trj,spe);
    for (j=0;j<numstep;++j)
      for (k=0;k<numdihed;++k) 
	speave[j*numdihed+k]=(i*speave[j*numdihed+k]+spe[j*numdihed+k])/(i+1);
  }
  fclose(inputfile1);

  /*****************************************************************************/
  /* if (flag=='v') {							       */
  /*   sumspe=(double *)gcemalloc(sizeof(double)*numdihed);		       */
  /*   for (i=0;i<numdihed;++i) sumspe[i]=0.0;				       */
  /*   for (i=0;(double)i/numstep/deltat/cv<1000.0;++i)			       */
  /*     for (j=0;j<numdihed;++j) sumspe[j]+=speave[i*numdihed+j];	       */
  /*   for (i=0;(double)i/numstep/deltat/cv<1000.0;++i)			       */
  /*     for (j=0;j<numdihed;++j) speave[i*numdihed+n]/=sumspe[j];	       */
  /* }									       */
  /* else {								       */
  /*   sumspe=(double *)gcemalloc(sizeof(double)*numdihed);		       */
  /*   for (i=0;i<numdihed;++i) sumspe[i]=0.0;				       */
  /*   for (i=0;(double)i/numstep/deltat/cv<1000.0;++i) {		       */
  /*     omg=(2.0*pi*i/numstep/deltat);					       */
  /*     for (j=0;j<numdihed;++j) {					       */
  /* 	sumspe[j]+=speave[i*numdihed+j]*omg*omg;			       */
  /*     }								       */
  /*   }								       */
  /*   for (i=0;(double)i/numstep/deltat/cv<1000.0;++i)			       */
  /*     for (j=0;j<numdihed;++j) 					       */
  /* 	speave[i*numdihed+j]=speave[i*numdihed+j]*omg*omg/sumspe[j];	       */
  /* }									       */
  /*****************************************************************************/

  outputfile=efopen(outputfilename,"w");
  pi=acos(-1.0);
  for (i=0;i<numstep;++i) {
    fprintf(outputfile,"%e ",(double)i/numstep/deltat/cv);
    for (j=0;j<numdihed;++j)
      if (flag=='v')
	fprintf(outputfile,"%e ",speave[i*numdihed+j]);
      else {
	omg=(2.0*pi*i/numstep/deltat);
	fprintf(outputfile,"%e ",speave[i*numdihed+j]*omg*omg);
      }
    fprintf(outputfile,"\n ");
  }
  fclose(outputfile);

  return 0;
}

void USAGE( char *progname) {
  printf("[-h] help\n");
  printf("[-c] coordinate mode\n" );
  printf("[-v] velocity mode\n");
  printf("[-n numdihed ]\n");
  printf("[-d deltat ]\n");
  printf("[-s numofsplit ]\n");
  printf("[-i numini ]\n");
  printf("[-f numfin ]\n");
  printf("[-t temp ]\n");
  printf("USAGE: %s inputfilename parmtopfilename outputfilename\n",progname);
}
