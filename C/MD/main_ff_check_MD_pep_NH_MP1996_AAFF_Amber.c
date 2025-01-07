#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "calcFFL_check.h"
#include "FFL.h"

#include "PTL.h"

#include "PTL.h"
#include "EF.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d;
  int numatom,numres;
  double dx=0.00001;
  int numspatom=10;

  double *crd,*mass,*vel;

  struct potential e;
  struct force f;
  double *pUMB,pUMB_t,**fUMB;

  int UMBdihedmode=OFF;

  double p_t=0.0,Etot;
  
  double pi;

  int *pairp;
  double *fcp,*dihe_equp;
  int nump;

  double f_e[3],f_LJ[3],f_e_14[3],f_LJ_14[3],f_d[3],f_a[3],f_b[3],f_UMB2[3];

  double x[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*parmfilename,*umbfilename;

  FILE *inputfile,*velfile,*parmfile,*umbfile;

  char *progname;

  int opt_idx=1;

  pi=acos(-1.0);

  struct option long_opt[] = {
    {"dx",1,NULL,'x'},
    {"nums",1,NULL,'s'},
    {"UMB",1,NULL,'u'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hx:s:u:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'x':
      dx=atof(optarg);
      break;
    case 's':
      numspatom=atoi(optarg);
      break;
    case 'u':
      umbfilename=optarg;
      UMBdihedmode=ON;
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
  inputfilename     = *argv;
  parmfilename      = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  j=0;
  for (i=0;i<numatom;++i) if (strncmp(AP.IGRAPH[i],"H",1)==0)  ++j;
  
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  if (UMBdihedmode==ON) {
    umbfile=efopen(umbfilename,"r");
    fscanf(umbfile,"%d",&nump);
    pairp=(int *)gcemalloc(sizeof(int)*nump*4);
    fcp=(double *)gcemalloc(sizeof(double)*nump);
    dihe_equp=(double *)gcemalloc(sizeof(double)*nump);
    pUMB=(int *)gcemalloc(sizeof(int)*nump);
    fUMB=(double **)gcemalloc(sizeof(double *)*numatom);
    for (i=0;i<numatom;++i)
      fUMB[i]=(double *)gcemalloc(sizeof(double)*3);
    for (i=0;i<nump;++i)
      fscanf(umbfile,"%d %d %d %d %lf %lf",&pairp[i*4+0],&pairp[i*4+1],&pairp[i*4+2],&pairp[i*4+3],&fcp[i],&dihe_equp[i]);
    fclose(umbfile);
    for (i=0;i<nump;++i) {
      for (j=0;j<4;++j) {
	pairp[i*4+j]-=1;
      }
      dihe_equp[i]=dihe_equp[i]*pi/180.0;
    }
  }

  ffL_set_calcffandforce(&e,&f);

  ffL_calcffandforce(crd,numatom,&e,&f);
  if (UMBdihedmode==ON) {
    pUMB_t=UMB_calc_dihetype_ff(crd,numatom,pairp,nump,fcp,dihe_equp,pUMB,fUMB);
  }

  ffL_calcffandforce_check(inputfilename,parmfilename,
			   numspatom,dx,
			   f_e,f_LJ,f_e_14,f_LJ_14,f_d,f_a,f_b);
  if (UMBdihedmode==ON) {
    UMB_calcffandforce_check(crd,numatom,pairp,nump,fcp,dihe_equp,
			     numspatom,dx,f_UMB2);
  }

  printf("f_es =( ");
  for (i=0;i<2;++i) {
    printf("%e ",f.f_e[numspatom*3+i]);
  }
  printf("%e ",f.f_e[numspatom*3+2]);
  printf(")\n");

  printf("f_es2=( ");
  for (i=0;i<2;++i) {
    printf("%e ",f_e[i]);
  }
  printf("%e ",f_e[2]);
  printf(")\n");

  printf("f_LJ =( ");
  for (i=0;i<2;++i) {
    printf("%e ",f.f_LJ[numspatom*3+i]);
  }
  printf("%e ",f.f_LJ[numspatom*3+2]);
  printf(")\n");

  printf("f_LJ2=( ");
  for (i=0;i<2;++i) {
    printf("%e ",f_LJ[i]);
  }
  printf("%e ",f_LJ[2]);
  printf(")\n");

  printf("f_14es =( ");
  for (i=0;i<2;++i) {
    printf("%e ",f.f_e_14[numspatom*3+i]);
  }
  printf("%e ",f.f_e_14[numspatom*3+2]);
  printf(")\n");

  printf("f_14es2=( ");
  for (i=0;i<2;++i) {
    printf("%e ",f_e_14[i]);
  }
  printf("%e ",f_e_14[2]);
  printf(")\n");

  printf("f_14LJ =( ");
  for (i=0;i<2;++i) {
    printf("%e ",f.f_LJ_14[numspatom*3+i]);
  }
  printf("%e ",f.f_LJ_14[numspatom*3+2]);
  printf(")\n");

  printf("f_14LJ2=( ");
  for (i=0;i<2;++i) {
    printf("%e ",f_LJ_14[i]);
  }
  printf("%e ",f_LJ_14[2]);
  printf(")\n");

  printf("f_d =( ");
  for (i=0;i<2;++i) {
    printf("%e ",f.f_d[numspatom*3+i]);
  }
  printf("%e ",f.f_d[numspatom*3+2]);
  printf(")\n");

  printf("f_d2=( ");
  for (i=0;i<2;++i) {
    printf("%e ",f_d[i]);
  }
  printf("%e ",f_d[2]);
  printf(")\n");

  printf("f_a =( ");
  for (i=0;i<2;++i) {
    printf("%e ",f.f_a[numspatom*3+i]);
  }
  printf("%e ",f.f_a[numspatom*3+2]);
  printf(")\n");

  printf("f_a2=( ");
  for (i=0;i<2;++i) {
    printf("%e ",f_a[i]);
  }
  printf("%e ",f_a[2]);
  printf(")\n");

  printf("f_b =( ");
  for (i=0;i<2;++i) {
    printf("%e ",f.f_b[numspatom*3+i]);
  }
  printf("%e ",f.f_b[numspatom*3+2]);
  printf(")\n");

  printf("f_b2=( ");
  for (i=0;i<2;++i) {
    printf("%e ",f_b[i]);
  }
  printf("%e ",f_b[2]);
  printf(")\n");

  if (UMBdihedmode==ON) {
    printf("f_UMB=( ");
    for (i=0;i<2;++i) {
      printf("%e ",fUMB[numspatom][i]);
    }
    printf("%e ",fUMB[numspatom][2]);
    printf(")\n");

    printf("f_UMB2=( ");
    for (i=0;i<2;++i) {
      printf("%e ",f_UMB2[i]);
    }
    printf("%e ",f_UMB2[2]);
    printf(")\n");
  }
    
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parmfilename\n",progname);
}
