#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PT.h"
#include "FF.h"
#include "EF.h"
#include "IOV2.h"
#include "EF.h"

#include "netcdf_mine.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,d;
  int numatom;
  double dx=0.0000001;

  double *crd,*crddx;
  double *ene,*f_check;
  struct potential e,e_dx;
  struct force f,f_dx;
  double fe[3];

  int flabb=ON,flaba=ON,flabd=ON,flabe=ON,flabl=ON;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*parmfilename;
  FILE *inputfile,*parmfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"b",0,NULL,'b'},
    {"a",0,NULL,'a'},
    {"d",0,NULL,'d'},
    {"e",0,NULL,'e'},
    {"l",0,NULL,'l'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"badelh",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'b':
      flabb=OFF;
      break;
    case 'a':
      flaba=OFF;
      break;
    case 'd':
      flabd=OFF;
      break;
    case 'e':
      flabe=OFF;
      break;
    case 'l':
      flabl=OFF;
      break;
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

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  parmfilename = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
 
  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddx=(double *)gcemalloc(sizeof(double)*numatom*3);
  f_check=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) f_check[i*3+j]=0.0;

  ff_set_calcffandforce(&e,&f);
  ff_set_calcffandforce(&e_dx,&f_dx);

  ff_calcffandforce(crd,numatom,&e,&f);

  if (flabb==ON) {
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k) {
 	memcpy(crddx,crd,sizeof(double)*numatom*3);
	crddx[j*3+k]+=dx;
	ff_calcffandforce(crddx,numatom,&e_dx,&f_dx);
	f_check[j*3+k]+=(e_dx.p_b_t-e.p_b_t)/dx*UNIT;
      }
    }
  }
  if (flaba==ON) {
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k) {
 	memcpy(crddx,crd,sizeof(double)*numatom*3);
	crddx[j*3+k]+=dx;
	ff_calcffandforce(crddx,numatom,&e_dx,&f_dx);
	f_check[j*3+k]+=(e_dx.p_a_t-e.p_a_t)/dx*UNIT;
      }
    }
  }
  if (flabd==ON) {
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k) {
 	memcpy(crddx,crd,sizeof(double)*numatom*3);
	crddx[j*3+k]+=dx;
	ff_calcffandforce(crddx,numatom,&e_dx,&f_dx);
	f_check[j*3+k]+=(e_dx.p_d_t-e.p_d_t)/dx*UNIT;
      }
    }
  }
  if (flabe==ON) {
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k) {
 	memcpy(crddx,crd,sizeof(double)*numatom*3);
	crddx[j*3+k]+=dx;
	ff_calcffandforce(crddx,numatom,&e_dx,&f_dx);
	f_check[j*3+k]+=(e_dx.p_e_t/*+1.0/1.2*e_dx.p_e_14_t*/-e.p_e_t/*-1.0/1.2*e.p_e_14_t*/)/dx*UNIT;
      }
    }
  }
  if (flabl==ON) {
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k) {
 	memcpy(crddx,crd,sizeof(double)*numatom*3);
	crddx[j*3+k]+=dx;
	ff_calcffandforce(crddx,numatom,&e_dx,&f_dx);
	f_check[j*3+k]+=(e_dx.p_LJ_t/*+0.5*e_dx.p_LJ_14_t*/-e.p_LJ_t/*-0.5*e.p_LJ_14_t*/)/dx*UNIT;
      }
    }
  }

  for (j=0;j<numatom;++j) {
    for (i=0;i<3;++i) fe[i]=0.0;
    if (flabb==ON) {
      for (k=0;k<3;++k) {
	fe[k]=f.f_b[j*3+k];
      }
    }
    if (flaba==ON) {
      for (k=0;k<3;++k) {
	fe[k]+=f.f_a[j*3+k];
      }
    }
    if (flabd==ON) {
      for (k=0;k<3;++k) {
	fe[k]+=f.f_d[j*3+k];
      }
    }
    if (flabe==ON) {
      for (k=0;k<3;++k) {
	fe[k]+=f.f_e[j*3+k]/*+1.0/1.2*f.f_e_14[j*3+k]*/;
      }
    }
    if (flabl==ON) {
      for (k=0;k<3;++k) {
	fe[k]+=f.f_LJ[j*3+k]/*+0.5*f.f_LJ_14[j*3+k]*/;
      }
    }
    printf("%2d: %8.3lf %8.3lf %8.3lf\n",j+1,fe[0],fe[1],fe[2]);
    printf("%2d: %8.3lf %8.3lf %8.3lf\n\n",j+1,f_check[j*3],f_check[j*3+1],f_check[j*3+2]);
  }
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("[-b] bond off \n");
  printf("[-a] angl off \n");
  printf("[-d] dihe off \n");
  printf("[-e] eles off \n");
  printf("[-l] LJ off \n");
  printf("%s [-h] inputfilename parmfilename \n",progname);
}
