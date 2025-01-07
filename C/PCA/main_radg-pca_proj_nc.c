
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
  int numatom,numatomp,numstep;
  double *traj,*cov,*eigenvalue;

  int MODE=AA,IOMODE=MD;

  char *inputfilename,*parmtopfilename,*vecterfilename;
  char *outputfilename;
  char *progname;

  FILE *inputfile,*parmtop,*vecterfile;
  FILE *outputfile;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  while((c=getopt(argc,argv,"Kbchs:"))!=-1) {
    switch(c) {
    case 'K':
      IOMODE=AMBER;
      break;
    case 'b':
      MODE=CA;
      break;
    case 'c':
      MODE=HV;
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  inputfilename      =  *argv;
  parmtopfilename    =  *++argv;
  vecterfilename     =  *++argv;
  outputfilename     =  *++argv;

  parmtop=efopen(parmtopfilename,"r");
  readParmtop(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  traj=mync_get_trj_aw(inputfilename,MODE,IOMODE,numatom,&numatomp,&numstep);
  //  traj=mync_get_trj_aw_b(inputfilename,inMODE,outMODE,IOMODE,numatom,&numatomp,&numstep);
  cov        = (double *)gcemalloc(sizeof(double)*numatomp*3*spdim);

  vecterfile=efopen(vecterfilename,"w");
  for (i=0;i<spdim;++i) {
    for (j=0;j<numatomp*3;++j) {
      fscanf(vecterfile,"%lf",&cov[i*numatomp*3+j]);
    }
  }
  fclose(vecterfile);

  RADG=RADG_calc_radg(crd,mass,numatomp);
  pca_proj_wdim(traj,cov,numstep,numatomp,spdim);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    for (j=0;j<spdim;++j) {
      fprintf(outputfile," %e",traj[i*numatomp*3+j]);
    }
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
