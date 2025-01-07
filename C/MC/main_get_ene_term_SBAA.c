
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
  int i,j,k,l;
  int flagT=ON,flagB=ON,flagA=ON,flagD=ON;
  int flagL=ON,flagL14=ON,flagE=ON,flagE14=ON,flagREST=ON;
  int flagmo=ON;
  int numini=1,numfin,interval=1;

  char *progname;
  char *inputfilename,*outfilename;
  FILE *outfile;

  double ene;
  struct my_netcdf_out_id_SBAAMCD nc_id_MCD;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"hi:f:v:"))!=-1) {
    switch(c) {
    case 'h':
      USAGE(progname);
      exit(1);
    case 'i':
      numini=atoi(optarg);
      break;
    case 'f':
      flagmo=OFF;
      numfin=atoi(optarg);
      break;
    case 'v':
      interval=atoi(optarg);
      break;
    default:
      USAGE(progname);
      exit(1);
    }
  }

  argc-=optind;
  argv+=optind;

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  outfilename  = *++argv;

  if (flagmo==ON)
    numfin=mync_get_present_step_SBAAMCD(inputfilename,&nc_id_MCD);

  if (numfin < numini){
    printf("error\n");
    exit(1);
  }

  if (numfin-numini < interval-1){
    printf("error\n");
    exit(1);
  }

  outfile=efopen(outfilename,"w");
  for (i=numini;i<numfin;i+=interval) {
    mync_open_inq_get_ene_SBAAMCD(inputfilename,i,1,i+1,&nc_id_MCD,&ene);
    fprintf(outfile,"%12.8e\n",ene);
  }
  fclose(outfile);

  return 0;
}

int USAGE(char *progname){
  printf("%s [-h] [-i numini] [-f numfin] [-v interval] inputfilename outfilename   \n",progname);
}
