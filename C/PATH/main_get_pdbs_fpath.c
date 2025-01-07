
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "PDB.h"
#include "PT.h"
#include "EF.h"

#include "netcdf_mine.h"

#define ON 1
#define OFF 0

#define AA 0
#define CA 1
#define HV 2

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,m;
  int MODE=AA;
  int numstep,numatom,interval=1;

  double *crd;

  int d;
  double f;

  PDBF PDB;

  double pi;
  double criteria;
  int *numpoints;
  double **path;
  
  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,outputfilename[100],*outputfilenamebase,*parmfilename,*progname;
  FILE *inputfile,*outputfile,*parmfile;

  while((c=getopt(argc,argv,"habci:"))!=-1) {
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
    case 'i':
      interval=atoi(optarg);
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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  parmfilename = *++argv;
  outputfilenamebase = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  PDB.numatom=AP.NATOM;
  PDB.PDBa=(PDBA *)gcemalloc(sizeof(PDBA)*AP.NATOM);
  readPDBdatafromParmtop(PDB);

  inputfile=efopen(inputfilename,"r");
  j=0;
  k=0;
  d = 1;
  while ( d != -1  )  {
    d=fscanf(inputfile,"%lf",&f);
    crd[j]=f;
    ++j;
    if (j==numatom*3) {
      if (k%interval==0) {
	for (l=0;l<numatom;++l) for (m=0;m<3;++m) PDB.PDBa[l].coord[m]=crd[l*3+m];
	sprintf(outputfilename,"%s_%d.pdb",outputfilenamebase,k);
	outputfile=efopen(outputfilename,"w");
	writPDB_wopt(outputfile,PDB,MODE);
	fclose(outputfile);
      }
      j=0;
      ++k;
    }
  }
  numstep=k;
  fclose(inputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-b] CA \n");
  printf("[-c] HV \n");
  printf("[-h] help \n");
  printf("%s [-b] [-c] [-h] inputfilename parmfilename outputfilenamebase\n",progname);
}

