#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "CHCOLLA.h"
#include "CB.h"
#include "CA.h"
#include "PT.h"
#include "EF.h"

#include "netcdf_mine.h"

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k;
  int numstep,numatom;
  double ave,var;

  int MODE=EXCH;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crd,*crdref;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;

  char *inputfilename,*outputfilename,*parmfilename,*progname;
  FILE *inputfile,*outputfile,*parmfile;

  while((c=getopt(argc,argv,"hH"))!=-1) {
    switch(c) {
    case 'H':
      MODE=INCH;
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
  parmfilename = *++argv;
  outputfilename = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdref=(double *)gcemalloc(sizeof(double)*numatom*3);
  mync_open_inq_get_sh_MCD(inputfilename,numatom,0,1,0+1,&nc_id_MCD,crd_nc);
  for (i=0;i<numatom;++i) 
    for (j=0;j<3;++j) crdref[i*3+j]=crd_nc[i][j];

  numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MCD);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
    for (j=0;j<numatom;++j) 
      for (k=0;k<3;++k) crd[j*3+k]=crd_nc[j][k];

    CHCO_diff_bond_lengh(crd,crdref,MODE,&ave,&var);
    fprintf(outputfile,"%d %e %e ",i,ave, var);
    CHCO_diff_bond_angle(crd,crdref,MODE,&ave,&var);
    fprintf(outputfile,"%e %e\n",ave, var);
  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-H] INCLUDE Hatom flag \n");
  printf("[-h] help \n");
  printf("%s [-H] [-h] inputfilename parmfilename outputfilename \n",progname);
}
