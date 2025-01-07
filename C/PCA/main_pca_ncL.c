
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PCA.h"
#include "PTL.h"
#include "EF.h"
#include "netcdf_mineL.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,s;
  int spdim=2;
  int numatom,numstep;
  double *traj,*cov,*eigenvalue;

  double *mass,massCA;

  int MODE=AA,INMODE=AA,NCMODE=MD;

  char *inputfilename,*parmtopfilename;
  char *outputfilename1,*outputfilename2,*outputfilename3;
  char *progname;

  FILE *inputfile, *parmtop;
  FILE *outputfile1,*outputfile2,*outputfile3;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  int opt_idx=1;

  struct option long_opt[] = {
    {"b",0,NULL,'b'},
    {"c",0,NULL,'c'},
    {"a",0,NULL,'a'},
    {"C",0,NULL,'C'},
    {"K",0,NULL,'K'},
    {"s",1,NULL,'s'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"bcaCKhs:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'K':
      NCMODE=AMBER;
      break;
    case 'C':
      INMODE=CA;
      break;
    case 'b':
      MODE=CA;
      break;
    case 'c':
      MODE=HV;
      break;
    case 'a':
      MODE=AA;
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

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  inputfilename      =  *argv;
  parmtopfilename    =  *++argv;
  outputfilename1    =  *++argv;
  outputfilename2    =  *++argv;
  outputfilename3    =  *++argv;

  parmtop=efopen(parmtopfilename,"r");
  readParmtopL(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  if (INMODE==CA) {
    j=0;
    for (i=0;i<numatom;++i) {
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
	++j; massCA=AP.AMASS[i];
      }
    }
    numatom=j;
  }

  mass=(double *)gcemalloc(sizeof(double)*numatom);
  if (INMODE==CA) for (i=0;i<numatom;++i) mass[i] = massCA;
  else for (i=0;i<numatom;++i) mass[i] = AP.AMASS[i];
  
  traj=myncL_get_trj_aw2(inputfilename,mass,NCMODE,numatom,&numstep);

  cov        = (double *)gcemalloc(sizeof(double)*numatom*3*numatom*3);
  eigenvalue = (double *)gcemalloc(sizeof(double)*numatom*3);

  pca_norm(traj,numstep,numatom);
  pca_covm(traj,numstep,numatom,cov);
  pca_diag(cov,eigenvalue,numatom);
  pca_proj(traj,cov,numstep,numatom);

  outputfile1=efopen(outputfilename1,"w");
  for (i=0;i<numstep;++i) {
    for (j=0;j<spdim;++j) fprintf(outputfile1," %e",traj[i*numatom*3+j]);
    fprintf(outputfile1,"\n");
  }
  fclose(outputfile1);

  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<spdim;++i) {
    for (j=0;j<numatom*3;++j) fprintf(outputfile2," %e",cov[i*numatom*3+j]);
    fprintf(outputfile2,"\n");
  }
  fclose(outputfile2);

  outputfile3=efopen(outputfilename3,"w");
  io_outputdata_f(outputfile3,numatom*3,eigenvalue);
  fclose(outputfile3);

  return 0;
}

void USAGE(char *progname) {
  printf("-b [CA]");
  printf("-c [exclude all H]");
  printf("-K [amber]");
  printf("%s inputfilename parmtopfilename outputfilename1(trj) outputfilename2(vec) outputfilename3(val)\n",progname);
}
