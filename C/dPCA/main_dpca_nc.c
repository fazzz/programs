
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PCA.h"
#include "dPCA.h"
#include "PTL.h"
#include "EF.h"
#include "IO.h"
#include "netcdf_mine.h"
//#include "netcdf_mineL.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,s;
  int spdim=2;
  int numatom,numdihed,numstep;
  double sum,sum2;
  double *dtrj,*sctrj,*cov,*eigenval,*dpca,tp;

  int MODE=AA,/*inMODE=AA,outMODE=CA,*/IOMODE=MD;

  char *inputfilename,*parmtopfilename;
  char *outputfilename,*outputfilename2,*outputfilename3;
  char *progname;

  FILE *inputfile, *parmtop;
  FILE *outputfile,*outputfile2,*outputfile3;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  int opt_idx=1;

  struct option long_opt[] = {
    {"K",0,NULL,'K'},
    {"s",1,NULL,'s'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"Khs:",long_opt,&opt_idx))!=-1) {
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

  if (argc < 7) {
    USAGE(progname);
    exit(1);
  }
  numdihed           =  atoi(*argv);
  numstep            =  atoi(*++argv);
  inputfilename      =  *++argv;
  parmtopfilename    =  *++argv;
  outputfilename     =  *++argv;
  outputfilename2    =  *++argv;
  outputfilename3    =  *++argv;
  
  parmtop=efopen(parmtopfilename,"r");
  readParmtopL(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  dtrj=(double *)gcemalloc(sizeof(double)*numstep*numdihed);
  sctrj=(double *)gcemalloc(sizeof(double)*numstep*numdihed*2);
  dpca=(double *)gcemalloc(sizeof(double)*numstep*numdihed*2);
  cov=(double *)gcemalloc(sizeof(double)*numdihed*2*numdihed*2);
  eigenval=(double *)gcemalloc(sizeof(double)*numdihed*2);

  for (i=0;i<numdihed*2;++i) for (j=0;j<numdihed*2;++j)  cov[i*numdihed*2+j]=0.0;
  inputfile=efopen(inputfilename,"r");
  //  io_scandtraj(inputfile,numstep,numdihed,dtrj);
  for (i=0;i<numstep;++i) {
    for (j=0;j<numdihed;++j) {
      fscanf(inputfile,"%lf",&dtrj[i*numdihed+j]);
    }
  }
  fclose(inputfile);

  dpca_norm(dtrj,sctrj,numstep,numdihed);
  dpca_covm(sctrj,numstep,numdihed,cov);
  dpca_diag(cov,eigenval,numdihed);
  dpca_proj(sctrj,dpca,cov,numstep,numdihed);
  
  sum=0.0;sum2=0;
  for (i=0;i<numdihed*2;++i) sum+=eigenval[i];
  outputfile=efopen(outputfilename,"w");
  fprintf(outputfile,"n eigenvalue(sigma^2) %% \n");
  for (i=0;i<numdihed*2;++i) {
    sum2+=eigenval[i];
    fprintf(outputfile,"%d %12.8lf %12.8lf \n",i,eigenval[i],sum2/sum*100.0);
  }
  fclose(outputfile);

  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<numdihed*2;++i) {
    for (j=0;j<spdim/*numdihed*2*/;++j)
      fprintf(outputfile2,"%lf ",cov[i*numdihed*2+j]);
    fprintf(outputfile2,"\n");
  }
  fclose(outputfile2);

  outputfile3=efopen(outputfilename3,"w");
  for (i=0;i<numstep;++i) {
    fprintf(outputfile3,"%d ",i);
    for (j=0;j<spdim/*numdihed*2*/;++j)
      fprintf(outputfile3,"%lf ",dpca[i*numdihed*2+j]);
    fprintf(outputfile3,"\n");
  }
  fclose(outputfile3);

  //  logfile=efopen("log_sc.txt","w");
  //  for (i=0;i<numstep;++i) {
  //    fprintf(logfile,"%d ",i);
  //    for (j=0;j<numdihed*2;++j)
  //      fprintf(logfile,"%lf ",sctrj[i*numdihed*2+j]);
  //    fprintf(logfile,"\n");
  //  }
  //  fclose(logfile);

  return 0;
}

void USAGE(char *progname) {
  printf("-b [CA]");
  printf("-c [exclude all H]");
  printf("-K [amber]");
  printf("%s inputfilename parmtopfilename outputfilename1(trj) outputfilename2(vec) outputfilename3(val)\n",progname);
}
