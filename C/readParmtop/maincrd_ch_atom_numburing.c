#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "PTL.h"
#include "EF.h"

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,numatom,d;
  double temp[3];
  double *crd;

  int NUM_CH_atoms;
  int *atom_numbering,*atom_numbering_mod;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *line;
  size_t len=0;

  char *parmfilename,*crdfilename,*crdfilemodname,*progname;
  FILE *parmfile,*crdfile,*crdfilemod;

  while((c=getopt(argc,argv,"h"))!=-1) {
    switch(c) {
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  parmfilename   = *argv;
  crdfilename    = *++argv;
  crdfilemodname = *++argv;
  NUM_CH_atoms   = atoi(*++argv);

  atom_numbering =(int *)gcemalloc(sizeof(int)*NUM_CH_atoms);
  atom_numbering_mod =(int *)gcemalloc(sizeof(int)*NUM_CH_atoms);

  for (i=0;i<NUM_CH_atoms;++i) {
    atom_numbering[i]=atoi(*++argv)-1;
    atom_numbering_mod[i]=atoi(*++argv)-1;
  }

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile);

  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  crdfile=efopen(crdfilename,"r");
  getline(&line,&len,crdfile);
  fscanf(crdfile,"%lf",&d);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      fscanf(crdfile,"%lf",&crd[i*3+j]);
    }
  }
  fclose(crdfile);

  for (i=0;i<NUM_CH_atoms;++i) {
    for (j=0;j<3;++j) {
      temp[j]=crd[atom_numbering[i]*3+j];
      crd[atom_numbering[i]*3+j]=crd[atom_numbering_mod[i]*3+j];
      crd[atom_numbering_mod[i]*3+j]=temp[j];
    }
  }
  
  crdfilemod=efopen(crdfilemodname,"w");
  fprintf(crdfilemod,"%s",line);
  fprintf(crdfilemod,"%d\n",numatom);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      if (crd[i*3+j]>=0.0)
	fprintf(crdfilemod,"  %-10.7lf",crd[i*3+j]);
      else
	fprintf(crdfilemod," -%-10.7lf",-crd[i*3+j]);
    }
    if ((i+1)%2 == 0)
      fprintf(crdfilemod,"\n");
  }
  fclose(crdfilemod);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] parmfilename crdfilename crdfilemodname NUM_CH_atoms numbering numbering_mod ...\n",progname);
}
