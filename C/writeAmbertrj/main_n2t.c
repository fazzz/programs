#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PDBL.h"

#include "TOPO.h"
#include "PT.h"
#include "EF.h"

#include "netcdf_mine.h"

#include "writeAmbertrj.h"

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,m;
  int numstep,numatom,interval=1;

  int flag;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crd;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id_MD;

  char *inputfilename,*outputfilename,*parmfilename,*progname,*proname;

  FILE *outputfile,*parmfile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"int",1,NULL,'i'},
    {"proname",1,NULL,'p'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  proname="    ";

  while((c=getopt_long(argc,argv,"hp:i:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'i':
      interval=atoi(optarg);
      break;
    case 'p':
      proname=optarg;
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
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MCD);
  //  makeAmbertrj(outputfilename,proname,outputfile);

  outputfile=efopen(outputfilename,"w");
  
  fprintf(outputfile," %4s\n",proname);

  for (i=0;i<numstep;++i) {
    mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);

    for (j=0;j<numatom;++j) for (k=0;k < 3; ++k) crd[j*3+k]=crd_nc[j][k];

    //    writeAmbertrj(crd,numatom,outputfile);
    l=1;
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k) {
	fprintf(outputfile," %7.3lf",crd[j*3+k]);
	if ((m=(l%10))==0) {
	  fprintf(outputfile,"\n");
	  flag=ON;
	}
	else 
	  flag=OFF;
	++l;
      }
    }
    if (flag==OFF) fprintf(outputfile,"\n");
  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[--int] interval \n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parmfilename outputfilename \n",progname);
}
