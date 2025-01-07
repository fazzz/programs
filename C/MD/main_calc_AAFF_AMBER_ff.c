#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "AAFF_AMBER_check.h"

#include "PTL.h"
#include "EF.h"
#include "NC.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numatom,numres;

  double f_e[3],f_LJ[3],f_e_14[3],f_LJ_14[3],f_d[3],f_a[3],f_b[3];
  double dx=0.00001;
  int numspatom=10,nums=3,numa=4;

  double *crd;

  struct potential e;
  struct force f;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*reffilename,*outputfilename,*parmfilename,*mapfilename;
  FILE *inputfile,*reffile,*outputfile,*parmfile,*mapfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"dx",1,NULL,'x'},
    {"nums",1,NULL,'s'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hx:s:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'x':
      dx=atof(optarg);
      break;
    case 's':
      numspatom=atoi(optarg);
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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  inputfilename     = *argv;
  parmfilename      = *++argv;
  outputfilename    = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  numres=AP.NRES;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  ffL_set_calcffandforce(&e,&f);
  ffL_calcffandforce(crd,numatom,&e,&f);

  AAFF_Amber_calcff_check(inputfilename,parmfilename,
			  numspatom,dx,
			  f_e,f_LJ,f_e_14,f_LJ_14,
			  f_d,f_a,f_b);

  printf("p_t=%e \n",e.p_t);
  printf("p_es =%e \n",e.p_e_t);
  printf("p_LJ =%e \n",e.p_LJ_t);
  printf("p_14e=%e \n",e.p_e_14_t);
  printf("p_14L=%e \n",e.p_LJ_14_t);
  printf("p_dih=%e \n",e.p_d_t);
  printf("p_ang=%e \n",e.p_a_t);
  printf("p_bon=%e \n",e.p_b_t);

  printf("f_es   =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f.f_e[numspatom*3+i]);
  }
  printf(")\n");

  printf("f_es   =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_e[i]);
  }
  printf(")\n");

  printf("f_LJ   =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f.f_LJ[numspatom*3+i]);
  }
  printf(")\n");

  printf("f_LJ   =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_LJ[i]);
  }
  printf(")\n");

  printf("f_14es =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f.f_e_14[numspatom*3+i]);
  }
  printf(")\n");

  printf("f_14es =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_e_14[i]);
  }
  printf(")\n");

  printf("f_14LJ =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f.f_LJ_14[numspatom*3+i]);
  }
  printf(")\n");

  printf("f_14LJ =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_LJ_14[i]);
  }
  printf(")\n");

  printf("f_d =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f.f_d[numspatom*3+i]);
  }
  printf(")\n");

  printf("f_d =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_d[i]);
  }
  printf(")\n");

  printf("f_a =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f.f_a[numspatom*3+i]);
  }
  printf(")\n");

  printf("f_a =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_a[i]);
  }
  printf(")\n");

  printf("f_b =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f.f_b[numspatom*3+i]);
  }
  printf(")\n");

  printf("f_b =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_b[i]);
  }
  printf(")\n");

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfilename parmfilename outputfilename\n",progname);
}


