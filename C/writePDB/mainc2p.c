#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <getopt.h>

#include "PDBL.h"

#include "TOPO.h"
#include "PT.h"
#include "EF.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numatom;

  int MODE=AA;

  PDBLF PDBL;

  char *line;
  size_t len=0;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crd;

  char *crdfilename,*pdbfilename,*parmfilename,*progname;
  FILE *crdfile,*pdbfile,*parmfile;

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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  crdfilename  = *argv;
  parmfilename = *++argv;
  pdbfilename  = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  PDBL.numatom=numatom;
  PDBL.PDBLa=(PDBLA *)gcemalloc(sizeof(PDBLA)*numatom);
  readPDBLdatafromParmtop(PDBL);

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  crdfile=efopen(crdfilename,"r");
  getline(&line,&len,crdfile);
  fscanf(crdfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(crdfile,"%lf",&crd[i*3+j]);
  fclose(crdfile);

  pdbfile=efopen(pdbfilename,"w");

  for (i=0;i<numatom;++i) for (j=0;j < 3; ++j) PDBL.PDBLa[i].coord[j]=crd[i*3+j];

  fprintf(pdbfile,"MODEL\n");
  writPDBL_wopt(pdbfile,PDBL,MODE);
  fprintf(pdbfile,"ENDMOD\n");
  fclose(pdbfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-H] (include Hydrogen) \n");
  printf("[--CA] (include Hydrogen) \n");
  printf("[-h] help \n");
  printf("%s [-h] crdfilename parmfilename pdbfilename \n",progname);
}
