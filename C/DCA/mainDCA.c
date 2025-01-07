#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "DCA.h"

#include "PT.h"
#include "FF.h"
#include "EF.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,d;
  int numstep,numatom,numclut,numtree;
  CLT *clt;
  AST *assemble_tree;
  double *Q,*frc;
  double *qacc,*qvel;

  double *crd,*mass;
  struct potential e;
  struct force f;

  double pi;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*outputfilename,*parmfilename,*clustfilename;
  FILE *inputfile,*outputfile,*parmfile,*clustfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"h",long_opt,&opt_idx))!=-1) {
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
  inputfilename   = *argv;
  clustfilename   = *++argv;
  parmfilename    = *++argv;
  outputfilename  = *++argv;

  pi=acos(-1.0);

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  clustfile=efopen(clustfilename,"r");
  clt=DCAp_clustscan(clustfile,&numclut);
  fclose(clustfile);

  DCAs_local_reference(clt,numclut,numatom,crd);
  DCAs_trans_Matrix(clt,numclut,numatom,crd);
  DCAs_inertia_matrix(clt,numclut,numatom,crd,mass);
  assemble_tree=DCAs_assemble_tree_testcase_2(numclut,&numtree);

  Q=(double *)gcemalloc(sizeof(double)*numclut);
  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  ff_set_calcffandforce(&e,&f);
  ff_calcffandforce(crd,numatom,&e,&f);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) frc[i*3+j]=f.f_e[i*3+j]+f.f_LJ[i*3+j]+f.f_e_14[i*3+j]+f.f_LJ_14[i*3+j];

  solverDCA(qacc,qvel,clt,Q,frc,assemble_tree,numclut,numatom,numtree);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfilename parmfilename  outputfilename\n",progname);
}
