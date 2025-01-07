
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PTL.h"

#include "netcdf_mine.h"

#define ON 0
#define OFF 1

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k;
  int numstep;

  double *crd;
  int numatom;

  int AMBERMODEflag=OFF;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id;

  char *line;
  size_t le=0;

  char *progname;
  char *inputfilename,*parmtopname;
  char *outputfilename;

  FILE *crdfile,*parmtop;
  FILE *outputfile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"Amber",0,NULL,'A'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hA",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'A':
      AMBERMODEflag=ON;
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
  inputfilename = *argv;
  parmtopname = *++argv;
  outputfilename = *++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtopL(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  if (AMBERMODEflag==ON) numstep=mync_get_present_step_AMBER(inputfilename,&nc_id);
  else numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MCD);

  outputfile=efopen(outputfilename,"w");
  fprintf(outputfile,"\n");
  for (i=0;i<numstep;++i) {
    if (AMBERMODEflag==ON) mync_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id,crd_nc);
    else mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
    for (j=0;j<numatom;++j)for (k=0;k<3;++k) crd[j*3+k]=crd_nc[j][k];

    for (j=0;j<numatom;++j) fprintf(outputfile,"%10.8lf %10.8lf %10.8lf \n",crd[j*3+0],crd[j*3+1],crd[j*3+2]);
  }
  fclose(outputfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("[--Amber] \n");
  printf("[--AA] \n");
  printf("%s [-h] inputfilename parmfilename outputfilename \n",progname);
}
