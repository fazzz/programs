#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "ABA.h"
#include "GOLMAA_hybrid_set.h"
#include "GOLMAA_hybrid.h"
#include "GOLMAA_hybrid_check.h"

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

  double f_natatt[3],f_repul[3];
  double dx=0.00001;
  int numspatom=52;

  double dt=0.001;
  double *crd,*refcrd,*mass;
  struct potential_GOLMAA_hybrid e;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*reffilename,*outputfilename,*parmfilename;
  FILE *inputfile,*reffile,*outputfile,*parmfile;

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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  inputfilename     = *argv;
  reffilename       = *++argv;
  parmfilename      = *++argv;
  outputfilename    = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  numres=AP.NRES;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);
  reffile=efopen(reffilename,"r");
  getline(&line,&len,reffile);
  fscanf(reffile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(reffile,"%lf",&refcrd[i*3+j]);
  fclose(reffile);

  GOLMAA_hybrid_ff_set_calcff(&e,refcrd,numatom,numres);

  GOLMAA_hyb_ff_calcff(crd,numatom,&e);
  GOLMAAff_hybrid_calcff_check(crd,numatom,&e,numspatom,dx,f_natatt,f_repul);

  printf("p_t=%e \n",e.p_t);
  printf("p_NC=%e \n",e.p_natatt_t);
  printf("p_NNC=%e \n",e.p_repul_t);

  printf("f_natt =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",e.f_natatt[numspatom][i]);
  }
  printf(")\n");

  printf("f_natt =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_natatt[i]);
  }
  printf(")\n");

  printf("f_repl =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",e.f_repul[numspatom][i]);
  }
  printf(")\n");

  printf("f_repl =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_repul[i]);
  }
  printf(")\n");

  outputfile=efopen(outputfilename,"w");
  fprintf(outputfile,"p_t=%e \n",e.p_t);
  fprintf(outputfile,"p_NC=%e \n",e.p_natatt_t);
  fprintf(outputfile,"p_NNC=%e \n",e.p_repul_t);
  fclose(outputfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfilename parmfilename outputfilename\n",progname);
}


