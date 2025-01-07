#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "PTL.h"
#include "FFL.h"
#include "netcdf_mineL.h"

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0;
  int numatom,numstep,interval=10;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crd;
  double crd_nc_in[MAXATOM][3],crd_nc_out[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD_in,nc_id_MCD_out;

  struct potential e;

  char *inputfilename,*outputfilename,*parmfilename;
  FILE *parmfile;

  char *progname;

  progname=*argv;

  while((c=getopt(argc,argv,"hi:"))!=-1) {
    switch(c) {
    case 'i':
      interval=atoi(optarg);
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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  parmfilename   = *++argv;
  outputfilename = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  numstep=myncL_get_present_step_MCD(inputfilename,&nc_id_MCD_in);
  myncL_create_def_MCD(outputfilename,numatom,&(nc_id_MCD_out));

  for (i=0;i<numstep;++i) {
    myncL_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD_in,crd_nc_in);

    if (i%interval==0) {
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc_out[j][k]=crd_nc_in[j][k];
      myncL_put_crd_ene_MCD(nc_id_MCD_out,l,crd_nc_out,e,0.0);
      ++l;
    }
  }
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[--int interval] set interval");
  printf("[-h] help \n");
  printf("%s [-i] [-h] inputfilename parmfilename outputfilename \n",progname);
}
