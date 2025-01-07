#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PDBL.h"
#include "writeAmberInput.h"

#include "PT.h"
#include "EF.h"

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0;
  int numatom;

  int MODE=AA;

  PDBLF PDBL;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crd;

  char *pdbfilename,*crdfilename,*parmfilename,*progname,*logfilename;
  FILE *pdbfile,*crdfile,*parmfile,*logfile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"CA",0,NULL,'C'},
    {"H",0,NULL,'H'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hCH",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'C':
      MODE=CA;
      break;
    case 'H':
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

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  pdbfilename = *argv;
  crdfilename = *++argv;

  pdbfile=efopen(pdbfilename,"r");
  readPDBLatomnum(pdbfile,&numatom);

  PDBL.PDBLa=(struct PDBLatom *)gcemalloc(sizeof(struct PDBLatom)*numatom);
  for (i=0;i<numatom;++i) {
    PDBL.PDBLa[i].HETEROflag=0;
    PDBL.PDBLa[i].serial=0;
    PDBL.PDBLa[i].ChainID=0;
    PDBL.PDBLa[i].resSeq=0;
    for (j=0;j<3;++j) PDBL.PDBLa[i].coord[j]=0.0;
    PDBL.PDBLa[i].occupancy=0.0;
    PDBL.PDBLa[i].tempfact=0.0;
  }

  pdbfile=efopen(pdbfilename,"r");
  readPDBL(pdbfile,PDBL,numatom);
  PDBL.numatom=numatom;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]=PDBL.PDBLa[i].coord[j];

  crdfile=efopen(crdfilename,"w");
  /**************************************************************/
  /* fprintf(crdfile,"\n");				        */
  /* fprintf(crdfile,"%d\n",numatom);			        */
  /* for (i=0;i<numatom;++i) {				        */
  /*   for (j=0;j<3;++j) fprintf(crdfile,"%8.4lf ",crd[i*3+j]); */
  /*   fprintf(crdfile,"\n");				        */
  /* }							        */
  /**************************************************************/
  writAmberInput(crdfile,numatom,crd,"ACE");
  fclose(crdfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-H] (include Hydrogen) \n");
  printf("[--CA] (include Hydrogen) \n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename outputfilename \n",progname);
}
