
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
  int i,j,k,l,m;
  int numatom;
  int flagmo=ON,flagl=OFF;
  int numini=1,numfin,interval=1,numl;

  char *progname;
  char *inputfilename,*outfilename;
  FILE *outfile;

  double crd[MAXATOM][3],*crd_data;
  struct my_netcdf_out_id_MCD nc_id_MCD;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"hi:f:v:l:"))!=-1) {
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
    case 'l':
      flagl=ON;
      numl=atoi(optarg);
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

  numatom=mync_get_numatom_MCD(inputfilename,&nc_id_MCD);
  crd_data=(double *)gcemalloc(sizeof(double)*numatom*3);

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
  if (flagl==OFF) {
    for (i=numini;i<numfin;i+=interval) {
      mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd);
      for (j=0;j<numatom;++j)
	for (k=0;k<3;++k)
	  crd_data[l]=crd[j][k];
      io_outputconf_Amberform(outfile,numatom,crd);
      fprintf(outfile,"\n");
    }
  }
  else {
    m=0;
    for (i=numfin-1;m<numl;i-=interval) {
      mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd);
      ++m;
      l=0;
      for (j=0;j<numatom;++j) {
	for (k=0;k<3;++k) {
	  crd_data[l]=crd[j][k];
	  ++l;
	}
      }
      io_outputconf_Amberform(outfile,numatom,crd);
      fprintf(outfile,"\n");
    }
  }
  fclose(outfile);

  return 0;
}

int USAGE(char *progname){
  printf("%s [-h] [-i numini] [-f numfin] [-v interval] [-l numla] inputfilename outfilename   \n",progname);
}
