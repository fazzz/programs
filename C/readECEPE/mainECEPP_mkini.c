#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ECEPE.h"
#include "PT.h"
#include "FF.h"

#include "EF.h"

#define ON  1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k;
  int numaom;
  int crdflag=OFF;

  double *co,*co_dummy,*dihed;

  struct ECEPE_parms ECEPE_p;
  struct pnb nb_p;

  char *progname;
  char *preofilename,*bd8filename,*coofilename,*outfilename;
  FILE *coofile,*outfile;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"hc"))!=-1) {
    switch(c) {
    case 'c':
      crdflag=ON;
      break;
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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  preofilename  = *argv;
  bd8filename = *++argv;
  coofilename = *++argv;
  outfilename = *++argv;

  read_ECEPE_parm(preofilename,bd8filename,&ECEPE_p,&nb_p);
  co=(double *)gcemalloc(sizeof(double)*ECEPE_p.NUMATM*3);
  dihed=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMVAR));
  co_dummy=(double *)gcemalloc(sizeof(double)*ECEPE_p.NUMATM*3);
  coofile=efopen(coofilename,"r");
  //  read_ECEPE_detail_coo(coofile,co,ECEPE_p.NUMATM);
  if (crdflag==OFF)
    read_ECEPE_coo(coofile,co_dummy,dihed,ECEPE_p.NUMATM);
  else {
    fscanf(coofile,"%d",&numaom);
    for (i=0;i<numaom;++i)
      for (j=0;j<3;++j)
	fscanf(coofile,"%lf",&co_dummy[i*3+j]);
  }    
  fclose(coofile);

  for (i=0;i<ECEPE_p.NUMATM;++i)
    for (j=0;j<3;++j) 
      co[(ECEPE_p.atom[i].katom-1)*3+j]=co_dummy[i*3+j];

  outfile=efopen(outfilename,"w");
  fprintf(outfile,"%d\n",ECEPE_p.NUMATM);
  for (i=0;i<ECEPE_p.NUMATM;++i) {
    for (j=0;j<3;++j) {
      fprintf(outfile,"%10.8lf ",co[i*3+j]);
    }
    fprintf(outfile,"\n");
  }
  fclose(outfile);

  return 0;
}

void USAGE(char *progname) {
  printf("-c -- crd flag ON\n");
  printf("-h -- help\n");
  printf("USAGE: %s profilename bd8filename coofilename(crdfilename) outfilename\n", progname);
}

