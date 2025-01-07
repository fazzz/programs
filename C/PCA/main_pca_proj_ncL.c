
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PCA.h"
#include "PTL.h"
#include "EF.h"
#include "IO.h"
#include "netcdf_mineL.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,s;
  int spdim=2;
  int numatom,numstep;
  double *traj,*cov,*eigenvalue;

  double *mass,massCA;

  int MODE=AA,inMODE=CA,IOMODE=MD;

  char *inputfilename,*parmtopfilename,*vecterfilename;
  char *outputfilename;
  char *progname;

  FILE *inputfile,*parmtop,*vecterfile;
  FILE *outputfile;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  int opt_idx=1;

  struct option long_opt[] = {
    {"b",0,NULL,'b'},
    {"c",0,NULL,'c'},
    {"K",0,NULL,'K'},
    {"C",0,NULL,'c'},
    {"s",1,NULL,'s'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"KCbchs:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'K':
      IOMODE=AMBER;
      break;
    case 'C':
      inMODE=CA;
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
  readParmtopL(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  if (inMODE==CA) {
    j=0; for (i=0;i<numatom;++i) { if (strncmp(AP.IGRAPH[i],"CA",2)==0){ ++j; massCA=AP.AMASS[i]; } }
    numatom=j;
  }

  mass=(double *)gcemalloc(sizeof(double)*numatom);
  if (inMODE==CA) for (i=0;i<numatom;++i) mass[i] = massCA;
  else for (i=0;i<numatom;++i) mass[i] = AP.AMASS[i];

  traj=myncL_get_trj_aw2(inputfilename,mass,IOMODE,numatom,&numstep);
  cov        = (double *)gcemalloc(sizeof(double)*numatom*3*spdim);

  vecterfile=efopen(vecterfilename,"w");
  for (i=0;i<spdim;++i) {
    for (j=0;j<numatom*3;++j) {
      fscanf(vecterfile,"%lf",&cov[i*numatom*3+j]);
    }
  }
  fclose(vecterfile);

  pca_proj_wdim(traj,cov,numstep,numatom,spdim);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    for (j=0;j<spdim;++j) {
      fprintf(outputfile," %e",traj[i*numatom*3+j]);
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
