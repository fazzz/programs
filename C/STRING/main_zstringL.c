
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <netcdf.h>
#include <math.h>

#include "STRING.h"
#include "CSI.h"
#include "PTL.h"
#include "FFL.h"
#include "EF.h"
#include "IO.h"

#include "netcdf_mineL.h"

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

  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }
  numiteration=atoi(*argv);
  numpoint=atoi(*++argv);
  inifilename  = *++argv;
  parmtopname  = *++argv;
  outputfilename=*++argv;
  outputfilename2=*++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtopL(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  path=(double *)gcemalloc(sizeof(double)*numatom*3*numpoint);
  inifile=efopen(inifilename,"r");
  for (i=0;i<sl;++i) getline(&line,&len,inifile);
  io_scantrj2(inifile,numatom,numpoint,path);
  fclose(inifile);

  path_evoluted=(double *)gcemalloc(sizeof(double)*numatom*3*numpoint);
  fe=(double *)gcemalloc(sizeof(double)*numatom*3*numpoint);
  pe=(double *)gcemalloc(sizeof(double)*numpoint);
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  ffL_set_calcffandforce(&e,&f);

  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<numiteration;++i) {
    for (j=0;j<numpoint*numatom*3;++j) fe[j]=0.0;
    for (j=0;j<numpoint;++j) pe[j]=0.0;

    for (j=0;j<numpoint;++j) {
      for (k=0;k<numatom*3;++k) crd[k]=path[j*numatom*3+k];
      ffL_calcffandforce(crd,numatom,&e,&f);
	if (flagB==ON) {
	  for (k=0;k<numatom*3;++k) fe[j*numatom*3+k]+=-f.f_b[k];
	  pe[j]+=e.p_b_t;
	}
	if (flagA==ON) {
	  for (k=0;k<numatom*3;++k) fe[j*numatom*3+k]+=f.f_a[k];
	  pe[j]+=e.p_a_t;
	}
	if (flagD==ON) {
	  for (k=0;k<numatom*3;++k) fe[j*numatom*3+k]+=f.f_d[k];
	  pe[j]+=e.p_d_t;
	}
	if (flagL==ON) {
	  for (k=0;k<numatom*3;++k)  fe[j*numatom*3+k]+=f.f_LJ[k]+f.f_LJ_14[k];
	  pe[j]+=e.p_LJ_t+e.p_LJ_14_t;
	}
	if (flagE==ON) {
	  for (k=0;k<numatom*3;++k)  fe[j*numatom*3+k]+=f.f_e[k]+f.f_e_14[k];
	  pe[j]+=e.p_e_t+e.p_e_14_t;
	}
    }
    z_string_Cartesian(path,path_evoluted,fe,numpoint,numatom,dt);
    if (i%outinterval==0) {
      io_outtrj2(outputfile,numatom,numpoint,path);
      for (j=0;j<numpoint;++j) fprintf(outputfile2,"%4d %12.8e \n",j,pe[j]);
    }
  }
  fclose(outputfile);
  fclose(outputfile2);
}

int USAGE(char *progname) {
  printf("USAGE:%s \n",progname);
  printf("[-B] BOND OFF [-A] ANGLE OFF [-D] DIHED OFF [-L] L-J [-E] elesta \n");
  printf("[-s sl] (dissmiss lines)    \n");
  printf("[-t dt]    \n");
  printf("[-o outinterval] \n");
  printf("[-h] help  \n");
  printf("%s  numiteration numpoint inifilename(path) parmtopname outputfilename(path) outputfilename2(ene)\n",progname);
}

 
