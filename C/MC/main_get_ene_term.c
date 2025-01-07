
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
  struct my_netcdf_out_id_MCD nc_id_MCD;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"hTBADLMECRi:f:v:"))!=-1) {
    switch(c) {
    case 'h':
      USAGE(progname);
      exit(1);
    case 'T':
      flagT=OFF;
      break;
    case 'B':
      flagB=OFF;
      break;
    case 'A':
      flagA=OFF;
      break;
    case 'D':
      flagD=OFF;
      break;
    case 'L':
      flagL=OFF;
      break;
    case 'M':
      flagL14=OFF;
      break;
    case 'E':
      flagE=OFF;
      break;
    case 'C':
      flagE14=OFF;
      break;
    case 'R':
      flagREST=OFF;
      break;
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
    numfin=mync_get_present_step_MCD(inputfilename,&nc_id_MCD);

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
    if (flagT==ON) {
      mync_open_inq_get_ene_MCD(inputfilename,i,1,i+1,0,&nc_id_MCD,&ene);
      fprintf(outfile,"%12.8e ",ene);
    }
    if (flagB==ON) {
      mync_open_inq_get_ene_MCD(inputfilename,i,1,i+1,1,&nc_id_MCD,&ene);
      fprintf(outfile,"%12.8e ",ene);
    }
    if (flagA==ON) {
      mync_open_inq_get_ene_MCD(inputfilename,i,1,i+1,2,&nc_id_MCD,&ene);
      fprintf(outfile,"%12.8e ",ene);
    }
    if (flagD==ON) {
      mync_open_inq_get_ene_MCD(inputfilename,i,1,i+1,3,&nc_id_MCD,&ene);
      fprintf(outfile,"%12.8e ",ene);
    }
    if (flagE==ON) {
      mync_open_inq_get_ene_MCD(inputfilename,i,1,i+1,4,&nc_id_MCD,&ene);
      fprintf(outfile,"%12.8e ",ene);
    }
    if (flagE14==ON) {
      mync_open_inq_get_ene_MCD(inputfilename,i,1,i+1,5,&nc_id_MCD,&ene);
      fprintf(outfile,"%12.8e ",ene);
    }
    if (flagL==ON) {
      mync_open_inq_get_ene_MCD(inputfilename,i,1,i+1,6,&nc_id_MCD,&ene);
      fprintf(outfile,"%12.8e ",ene);
    }
    if (flagL14==ON) {
      mync_open_inq_get_ene_MCD(inputfilename,i,1,i+1,7,&nc_id_MCD,&ene);
      fprintf(outfile,"%12.8e ",ene);
    }
    if (flagREST==ON) {
      mync_open_inq_get_ene_MCD(inputfilename,i,1,i+1,8,&nc_id_MCD,&ene);
      fprintf(outfile,"%12.8e ",ene);
    }
    fprintf(outfile,"\n",ene);
  }
  fclose(outfile);

  return 0;
}

int USAGE(char *progname){
  printf("%s [-h] [-T] [-B] [-A] [-D] [-L] [-M] [-E] [-C] [-R] [-i numini] [-f numfin] [-v interval] inputfilename outfilename   \n",progname);
}
