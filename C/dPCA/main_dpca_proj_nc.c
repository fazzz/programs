
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>

#include "PCA.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"
#include "netcdf_mine.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,s;
  int spdim=2;
  int numatom,numdihed,numstep;
  double *dtrj,*sctrj,*cov,*eigenvalue,*dpca;

  int MODE=AA,IOMODE=MD;

  char *inputfilename,*parmtopfilename,*vecterfilename;
  char *outputfilename;
  char *progname;

  FILE *inputfile,*parmtop,*vecterfile;
  FILE *outputfile;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  while((c=getopt(argc,argv,"Khs:"))!=-1) {
    switch(c) {
    case 'K':
      IOMODE=AMBER;
      break;
    case 's':
      spdim=atoi(optarg);
      break;
    case 'h':
      USAGE(progname);
      break;
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=argv[0];
  argc-=optind;
  argv+=optind;

  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }
  numdihed           =  atoi(*argv);
  numstep            =  atoi(*++argv);
  inputfilename      =  *++argv;
  parmtopfilename    =  *++argv;
  vecterfilename     =  *++argv;
  outputfilename     =  *++argv;

  parmtop=efopen(parmtopfilename,"r");
  readParmtop(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  dtrj=(double *)gcemalloc(sizeof(double)*numstep*numdihed);
  sctrj=(double *)gcemalloc(sizeof(double)*numstep*numdihed*2);
  //  dpca=(double *)gcemalloc(sizeof(double)*numstep*numdihed*2);
  cov=(double *)gcemalloc(sizeof(double)*numdihed*2*spdim);
  // eigenval=(double *)gcemalloc(sizeof(double)*numdihed*2);
  dpca=(double *)gcemalloc(sizeof(double)*numstep*spdim);

  inputfile=efopen(inputfilename,"r");
  for (i=0;i<numstep;++i) {
    for (j=0;j<numdihed;++j) {
      fscanf(inputfile,"%lf",&dtrj[i*numdihed+j]);
    }
  }
  fclose(inputfile);

  for (i=0;i<numdihed*2*spdim;++i) cov[i]=0.0;
  vecterfile=efopen(vecterfilename,"r");
  for (i=0;i<numdihed*2;++i) {
    for (j=0;j<spdim;++j) {
      fscanf(vecterfile,"%lf",&cov[i*spdim+j]);
    }
  }
  fclose(vecterfile);

  dpca_norm(dtrj,sctrj,numstep,numdihed);
  dpca_proj_wdim(sctrj,dpca,cov,numstep,numdihed,spdim);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    //    fprintf(outputfile,"%d ",i);
    for (j=0;j<spdim;++j)  fprintf(outputfile," %e",dpca[i*spdim+j]);
    fprintf(outputfile,"\n");
  }
  fclose(outputfile);

  return 0;
}

void USAGE(char *progname) {
  printf("-b [CA]");
  printf("-c [exclude all H]");
  printf("-K [amber]");
  printf("-s [numdim]");
  printf("-h [help]");
  printf("%s inputfilename parmtopfilename vecterfilename outputfilename\n",progname);
}
