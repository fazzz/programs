#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "netcdf_mine.h"

#include "PDB.h"
#include "MB.h"
#include "PTL.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0

int usage(void);

int main(int argc, char *argv[]) {
  int i,j;
  int *numdihed,nt=1,nd;
  int flag='P';
  int **adpairs;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *parmfilename;
  char *progname;
  FILE *parmfile;

  while((c=getopt(argc,argv,"hPOK"))!=-1) {
    switch(c) {
    case 'P':
      flag='P';
      nt=1;
      break;
    case 'O':
      flag='O';
      nt=2;
      break;
    case 'K':
      flag='K';
      nt=5;
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

  if (argc < 1) {
    USAGE(progname);
    exit(1);
  }
  parmfilename  =  *argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 

  adpairs=(int **)gcemalloc(sizeof(int *)*5);
  adpairs[0]=(int *)gcemalloc(sizeof(int)*4); // PHI,PSI
  adpairs[1]=(int *)gcemalloc(sizeof(int)*4); // OMEGA
  adpairs[2]=(int *)gcemalloc(sizeof(int)*4); // KI
  adpairs[3]=(int *)gcemalloc(sizeof(int)*8); // ACE
  adpairs[4]=(int *)gcemalloc(sizeof(int)*8); // NME
  numdihed=(int *)gcemalloc(sizeof(int)*5);  
  readdihedpairsL(adpairs,numdihed);
  
  nd=0;
  for (i=0;i<nt;++i) {
    nd+=numdihed[i];
  }

  printf("%d \n",nd);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-P] phi psi \n");
  printf("[-O] phi psi omega \n");
  printf("[-K] phi psi omega kai \n");
  printf("%s inputfilename parmfilename outputfilename \n",progname);
}

 
