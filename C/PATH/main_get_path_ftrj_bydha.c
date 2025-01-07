
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "PATH.h"
#include "rmsd.h"
#include "PT.h"
#include "EF.h"

#include "netcdf_mine.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int MODE=AA;
  int numpath,numatom,interval=1;

  double *crdref1,*crdref2;

  double *dha_spe;
  int *atomp;
  int numdha;

  int d;
  double f;

  double pi;
  double criteria;
  int *numpoints;
  double **path;
  
  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,outputfilename,*outputfilenamebase,*parmfilename,*reffilename1,*reffilename2,*condfilename,*progname;
  FILE *outputfile,*parmfile,*reffile1,*reffile2,*condfile;

  while((c=getopt(argc,argv,"habc"))!=-1) {
    switch(c) {
    case 'a':
      MODE=AA;
      break;
    case 'b':
      MODE=CA;
      break;
    case 'c':
      MODE=HV;
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

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  criteria = atof(*argv)/180*pi;
  inputfilename = *++argv;
  condfilename = *++argv;
  parmfilename = *++argv;
  outputfilenamebase = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  condfile=efopen(condfilename,"r");
  fscanf(condfile,"%d",&numdha);
  atomp=(int *)gcemalloc(sizeof(int)*numdha*4);
  dha_spe=(double *)gcemalloc(sizeof(double)*numdha*2);
  for (i=0;i<numdha;++i) {
    for (j=0;j<4;++j) {
      fscanf(condfile,"%d",&d);
      atomp[i*4+j]=d-1;
    }
    for (j=0;j<2;++j) {
      fscanf(condfile,"%lf",&f);
      dha_spe[i*2+j]=f/180.0*pi;
      if (dha_spe[j] > pi) dha_spe[j]-=2.0*pi;
      else if (dha_spe[j] <  -1.0*pi) dha_spe[j]+=2.0*pi;
    }
  }
  fclose(condfile);

  path=(double **)gcemalloc(sizeof(double *)*1);
  path[0]=(double *)gcemalloc(sizeof(double)*numatom*3);
  numpoints=(int *)gcemalloc(sizeof(int)*1);
  numpath=get_path_ftrj_bydha(inputfilename,outputfilenamebase,crdref1,crdref2,criteria,dha_spe,atomp,numdha,path,numpoints,MODE,interval);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-b] CA \n");
  printf("[-c] HV \n");
  printf("[-h] help \n");
  printf("%s [-b] [-c] [-h] criteria inputfilename condfilename parmfilename outputfilenamebase\n",progname);
}

