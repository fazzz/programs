
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <netcdf.h>

#include "EF.h"
#include "netcdf_mine.h"

#define ON 1
#define OFF 0

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int nstep;

  char *progname;
  char *inputfilename;

  struct my_netcdf_out_id_MCD nc_id_MCD;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"h"))!=-1) {
    switch(c) {
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

  if (argc < 1) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;

  nstep=mync_get_present_step_MCD(inputfilename,&nc_id_MCD);

  printf("%d\n",nstep);

  return 0;
}

int USAGE(char *progname){
  printf("%s inputfilename\n",progname);
}
