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

#include "netcdf_mineL.h"

#include "writeAmbertrj.h"

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,m;
  int numstep,numatom,interval=1;

  int AMBncflag=OFF,nameflag=ON;
  int flag,flag_spe_numatom=OFF;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crd;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id;

  char *inputfilename,*outputfilename,*parmfilename,*progname,*proname;

  FILE *outputfile,*parmfile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"int",1,NULL,'i'},
    {"proname",1,NULL,'p'},
    {"numatom",1,NULL,'a'},
    {"nameflagOFF",0,NULL,'n'},
    {"Amber",0,NULL,'K'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  proname="    ";

  while((c=getopt_long(argc,argv,"hnKp:i:a:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'i':
      interval=atoi(optarg);
      break;
    case 'p':
      proname=optarg;
      break;
    case 'a':
      numatom=atoi(optarg);
      flag_spe_numatom=ON;
      break;
    case 'n':
      nameflag=OFF;
      break;
    case 'K':
      AMBncflag=ON;
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
  if (flag_spe_numatom==OFF)
    numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  if (AMBncflag==ON) numstep=myncL_get_present_step_AMBER(inputfilename,&nc_id);
  else numstep=myncL_get_present_step_MCD(inputfilename,&nc_id_MCD);

  outputfile=efopen(outputfilename,"w");
  
  if (nameflag==ON) {
    fprintf(outputfile," %4s\n",proname);
  }

  for (i=0;i<numstep;++i) {
    if (AMBncflag==ON) myncL_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id,crd_nc);
    else myncL_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);

    for (j=0;j<numatom;++j) for (k=0;k < 3; ++k) crd[j*3+k]=crd_nc[j][k];

    //    writeAmbertrj(crd,numatom,outputfile);
    l=1;
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k) {
	//	fprintf(outputfile," %7.3lf",crd[j*3+k]);
	fprintf(outputfile," %12.8lf",crd[j*3+k]);
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
