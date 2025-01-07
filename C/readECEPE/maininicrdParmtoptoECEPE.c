#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ECEPE.h"
#include "PT.h"

#include "EF.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,m,n;
  int numatom,numres;

  int *indexkatom,*indexjatom;

  double *crd;

  struct ECEPE_parms ECEPE_p;
  struct pnb nb_p;

  char *progname;
  char *preofilename,*bd8filename,*parmtopfilename,*crdfilename;
  FILE *parmtopfile,*crdfile;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
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

  argc-=optind;
  argv+=optind;

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  preofilename    = *argv;
  bd8filename     = *++argv;
  parmtopfilename = *++argv;
  crdfilename     = *++argv;
  //  outputfilename = *++argv;

  read_ECEPE_parm(preofilename,bd8filename,&ECEPE_p,&nb_p);

  parmtopfile=efopen(parmtopfilename,"r");
  readParmtop(parmtopfile);
  fclose(parmtopfile);

  indexkatom=(int *)gcemalloc(sizeof(int)*AP.NATOM);
  indexjatom=(int *)gcemalloc(sizeof(int)*AP.NATOM);

  for (i=0;i<AP.NATOM;++i) {
    indexkatom[i]=-1;
    indexjatom[i]=-1;
  }

  numres=-1;
  for (i=0;i<AP.NATOM;++i) {
    if (i>=AP.IPRES[numres+1]-1) ++numres;
    for (j=0;j<ECEPE_p.NUMATM;++j) {
      if (strncmp(AP.IGRAPH[i],ECEPE_p.atom[j].name_atom,2)==0 
	  && strncmp(AP.LABERES[numres],ECEPE_p.atom[j].name_res,3)==0 ) {
	indexkatom[i]=ECEPE_p.atom[j].katom;
	indexjatom[i]=ECEPE_p.atom[j].jatom;
	break;
      }
    }
  }

  return 0;
}

void USAGE(char *progname) {
  printf("-h -- help\n");
  printf("USAGE: %s profilename bd8filename crdfilename outputfilenam\n", progname);
}

