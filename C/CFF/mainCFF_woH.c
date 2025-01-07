
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "STRING.h"
#include "CSI.h"
#include "PT.h"
#include "FF.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int flag=ON,flagB=ON,flagA=ON,flagD=ON,flagL=ON,flagE=ON;

  int sl=0;

  int outinterval=1;

  int numpoint,numiteration;
  int numatom,numpara,numnb,num14;


  double dt=0.01;

  char *line;
  size_t len=0;

  char *progname;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *path,*path_evoluted;
  double *fe,*pe;

  struct potential e;
  struct force f;
  double *crd;

  char *inifilename,*parmtopname,*outputfilename,*outputfilename2;
  FILE *inifile,*parmtop,*outputfile,*outputfile2;

  progname=argv[0];

  while((c=getopt(argc,argv,"BADLEhs:t:o:"))!=-1) {
    switch(c) {
    case 'B':
      flagB=OFF;
      break;
    case 'A':
      flagA=OFF;
      break;
    case 'D':
      flagD=OFF;
      break;
    case 'L':
      flagL=OFF;
      break;
    case 'E':
      flagE=OFF;
      break;
    case 's':
      sl=atoi(optarg);
      break;
    case 't':
      dt=atof(optarg);
      break;
    case 'o':
      outinterval=atoi(optarg);
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  numpoint=atoi(*argv);
  inifilename  = *++argv;
  parmtopname  = *++argv;
  outputfilename=*++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtop(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  path=(double *)gcemalloc(sizeof(double)*numatom*3*numpoint);
  inifile=efopen(inifilename,"r");
  for (i=0;i<sl;++i) getline(&line,&len,inifile);
  io_scantrj2(inifile,numatom,numpoint,path);
  fclose(inifile);

  fe=(double *)gcemalloc(sizeof(double)*numatom*3*numpoint);
  pe=(double *)gcemalloc(sizeof(double)*numpoint);
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  ff_set_calcffandforce(&e,&f);

  outputfile=efopen(outputfilename,"w");
  for (j=0;j<numpoint*numatom*3;++j) fe[j]=0.0;
  for (j=0;j<numpoint;++j) pe[j]=0.0;

  for (j=0;j<numpoint;++j) {
    for (k=0;k<numatom*3;++k) crd[k]=path[j*numatom*3+k];
    ff_calcffandforce_woH(crd,numatom,&e,&f);
    fprintf(outputfile,"%4d %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",j,e.p_e_t+e.p_LJ_t+e.p_e_14_t+e.p_LJ_14_t+e.p_d_t+e.p_a_t+e.p_b_t,e.p_e_t,e.p_LJ_t,e.p_e_14_t,e.p_LJ_14_t,e.p_d_t,e.p_a_t,e.p_b_t);
  }
  fclose(outputfile);
}

int USAGE(char *progname) {
  printf("USAGE:%s \n",progname);
  printf("[-B] BOND OFF [-A] ANGLE OFF [-D] DIHED OFF [-L] L-J [-E] elesta \n");
  printf("[-s sl] (dissmiss lines)    \n");
  printf("[-t dt]    \n");
  printf("[-o outinterval] \n");
  printf("[-h] help  \n");
  printf("%s  numpoint inifilename(path) parmtopname  outputfilename(ene)\n",progname);
}

 
