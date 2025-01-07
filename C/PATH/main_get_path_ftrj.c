
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

  double criteria;
  int *numpoints;
  double **path;
  
  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,outputfilename,*outputfilenamebase,*parmfilename,*reffilename1,*reffilename2,*progname;
  FILE *outputfile,*parmfile,*reffile1,*reffile2;

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

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }
  criteria = atof(*argv);
  inputfilename = *++argv;
  reffilename1 = *++argv;
  reffilename2 = *++argv;
  parmfilename = *++argv;
  outputfilenamebase = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  crdref1=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdref2=(double *)gcemalloc(sizeof(double)*numatom*3);

  reffile1=efopen(reffilename1,"r");
  io_scanconf_Amber_ini(reffile1,numatom,crdref1);
  fclose(reffile1);

  reffile2=efopen(reffilename2,"r");
  io_scanconf_Amber_ini(reffile2,numatom,crdref2);
  fclose(reffile2);

  path=(double **)gcemalloc(sizeof(double *)*1);
  path[0]=(double *)gcemalloc(sizeof(double)*numatom*3);
  numpoints=(int *)gcemalloc(sizeof(int)*1);
  numpath=get_path_ftrj(inputfilename,outputfilenamebase,crdref1,crdref2,criteria,path,numpoints,MODE,interval);

  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-b] CA \n");
  printf("[-c] HV \n");
  printf("[-h] help \n");
  printf("%s [-b] [-c] [-h]   criteria inputfilename reffile1name reffile2name parmfilename outputfilenamebase\n",progname);
}

