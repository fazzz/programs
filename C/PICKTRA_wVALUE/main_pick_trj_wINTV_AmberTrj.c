#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "TOPO.h"
#include "LA.h"
#include "mymath.h"

#include "MB.h"
#include "PTL.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0

int usage(void);

int main(int argc, char *argv[]) {
  int i,j,k,l,d,dummy;
  int numatom,numstep=100;

  double *crd;

  int TLflag=ON;
  int SetNUMCRDflag=OFF;

  int numcrd=0,NUMCRDtotall=10;
  int initialstep=0,interval=1;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*parmfilename;
  char *crdfilenamebase,crdfilename[100];
  char *progname;
  FILE *inputfile,*parmfile,*crdfile;

  while((c=getopt(argc,argv,"hTs:v:n:i:"))!=-1) {
    switch(c) {
    case 'T':
      TLflag=OFF;
      break;
    case 's':
      numstep=atoi(optarg);
      break;
    case 'v':
      interval=atoi(optarg);
      break;
    case 'n':
      NUMCRDtotall=atoi(optarg);
      SetNUMCRDflag=ON;
      break;
    case 'i':
      initialstep=atoi(optarg);
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
  inputfilename   = *argv;
  parmfilename    = *++argv; 
  crdfilenamebase = *++argv;
  
  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;
  
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfile=efopen(inputfilename,"r");
  if (TLflag==ON)   getline(&line,&len,inputfile);

  if (initialstep>0)
    for (i=0;i<initialstep;++i) for (j=0;j<numatom;++j) for (k=0;k<3;++k) fscanf(inputfile,"%lf",&crd[j*3+k]);

  if (SetNUMCRDflag==ON) interval=(int)(numstep/NUMCRDtotall);

  if (interval==0) {
    printf("ERORR interval must be integer:");
    exit(0);
  }

  for (i=initialstep;i<numstep;++i) {
    for (j=0;j<numatom;++j) for (k=0;k<3;++k) fscanf(inputfile,"%lf",&crd[j*3+k]);

    if ( (i%interval) == 0 || i==initialstep ) {
      ++numcrd;
      sprintf(crdfilename,"%s_%d",crdfilenamebase,numcrd);
      crdfile=efopen(crdfilename,"w");
      
      io_outputconf_Amberform(crdfile,numatom,crd);
      fclose(crdfile);
    }
  }

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-T ] T(op)L(ine)flag OFF \n");
  printf("[-s numstep] specify number of steps \n");
  printf("[-v interval] specify interval steps \n");
  printf("[-i initial] specify initial steps \n");
  printf("[-h ] help \n");
  printf("%s inputfilename parmfilename crdfilenamebase \n",progname);
}
