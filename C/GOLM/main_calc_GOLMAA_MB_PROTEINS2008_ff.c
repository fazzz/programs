#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"
#include "GOLMAA_MB_PROTEINS2008.h"
#include "GOLMAA_MB_PROTEINS2008_check.h"

#include "PTL.h"
#include "EF.h"
#include "NC.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,a;
  int numatom,numres;

  double f_MB[3],f_e1[3],f_e2[3];
  double dx=0.00001;
  int numspatom=10,nums=3,numa=4;

  double *crd,*refcrd1,*refcrd2;
  double d=1.0,de=1.0,d2;

  double ep=ep_natatt_PROTEINS2008;

  int nibnum=3,criteria=criteria_PROTEINS2008;

  struct potential e;
  struct force f;
  struct potential_GOLMAA_MB_PROTEINS2008 e_GOLM;

  double x[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*reffilename1,*reffilename2,*outputfilename,*parmfilename,*mapfilename;
  FILE *inputfile,*reffile1,*reffile2,*outputfile,*parmfile,*mapfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"dx",1,NULL,'x'},
    {"nums",1,NULL,'s'},
    {"numa",1,NULL,'a'},
    {"ep",1,NULL,'e'},
    {"cutoff",1,NULL,'c'},
    {"de",1,NULL,'d'},
    {"d",1,NULL,'2'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hx:s:a:e:c:d:2:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'x':
      dx=atof(optarg);
      break;
    case 's':
      numspatom=atoi(optarg);
      break;
    case 'a':
      numa=atoi(optarg);
      break;
    case 'e':
      ep=atof(optarg);
      break;
    case 'c':
      criteria=atof(optarg);
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  inputfilename     = *argv;
  reffilename1      = *++argv;
  reffilename2      = *++argv;
  parmfilename      = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  numres=AP.NRES;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd1=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd2=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&a);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  reffile1=efopen(reffilename1,"r");
  getline(&line,&len,reffile1);
  fscanf(reffile1,"%d",&a);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(reffile1,"%lf",&refcrd1[i*3+j]);
  fclose(reffile1);

  reffile2=efopen(reffilename2,"r");
  getline(&line,&len,reffile2);
  fscanf(reffile2,"%d",&a);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(reffile2,"%lf",&refcrd2[i*3+j]);
  fclose(reffile2);

  ffL_set_calcffandforce(&e,&f);
  GOLMAA_MB_PROTEINS2008_ff_calcff_set(&e_GOLM,refcrd1,refcrd2,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);

  d2=d*d;
  GOLMAA_MB_PROTEINS2008_ff_calcff(crd,numatom,de,d2,&e_GOLM);
  GOLMAA_MB_PROTEINS2008_calcff_check(inputfilename,reffilename1,reffilename2,parmfilename,
				      numspatom,dx,d,de,ep,nibnum,criteria,f_MB,f_e1,f_e2);

  printf("p_t=%e \n",e_GOLM.p_MB);

  printf("f_MB =( ");
  for (i=0;i<2;++i) {
    printf("%e ",e_GOLM.f_MB[numspatom][i]);
  }
  printf("%e ",e_GOLM.f_MB[numspatom][2]);
  printf(")\n");

  printf("f_MB2=( ");
  for (i=0;i<2;++i) {
    printf("%e ",f_MB[i]);
  }
  printf("%e ",f_MB[2]);
  printf(")\n");

  printf("f_e1 =( ");
  for (i=0;i<2;++i) {
    printf("%e ",e_GOLM.e1.f_t[numspatom][i]);
  }
  printf("%e ",e_GOLM.e1.f_t[numspatom][2]);
  printf(")\n");

  printf("f_e12=( ");
  for (i=0;i<2;++i) {
    printf("%e ",f_e1[i]);
  }
  printf("%e ",f_e1[2]);
  printf(")\n");

  printf("f_e2 =( ");
  for (i=0;i<2;++i) {
    printf("%e ",e_GOLM.e2.f_t[numspatom][i]);
  }
  printf("%e ",e_GOLM.e2.f_t[numspatom][2]);
  printf(")\n");

  printf("f_e22=( ");
  for (i=0;i<2;++i) {
    printf("%e ",f_e2[i]);
  }
  printf("%e ",f_e2[2]);
  printf(")\n");

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfilename parmfilename outputfilename\n",progname);
}

