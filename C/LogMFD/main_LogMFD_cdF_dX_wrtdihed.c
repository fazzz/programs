
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "netcdf_mineL.h"
#include "PTL.h"
#include "FFLc.h"
#include "EF.h"

#include "LogMFD_cMF.h"

#define ON 0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k;
  int na1,na2;
  int numatom,numstep;

  int flag=OFF;

  double dF_dX,dF_dX_a;

  double crd_nc[MAXATOM][3],*crd;
  struct my_netcdf_out_id_AMBER nc_id_MD;

  double pi;

  char *inputfilename,*outputfilename,*parmfilename;
  FILE *outputfile,*parmfile;
  
  char *line;
  size_t len=0;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;

  struct AmberParmL ap;

  struct potential e;
  struct force f;

  pi=acos(-1.0);
  
  struct option long_opt[] = {
    {"O",0,NULL,'o'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"ho",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'o':
      flag=ON;
      break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);  exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  outputfilename = *++argv;
  parmfilename   = *++argv;
  na1 = atoi(*++argv)-1;
  na2 = atoi(*++argv)-1;

  numstep=myncL_get_present_step_AMBER(inputfilename,&nc_id_MD);

  parmfile=efopen(parmfilename,"r");
  readParmtopLb(parmfile,&ap);
  fclose(parmfile);

  numatom=ap.NATOM;
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  ffLc_set_calcffandforce(&e,&f,ap);

  for (i=0;i<numstep;++i) {
    myncL_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id_MD,crd_nc);

    for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd[j*3+k]=crd_nc[j][k];

    dF_dX=LogMFD_cdF_dX_wrtdihed(crd,&e,&f,ap,numatom,na1,na2);

    dF_dX_a=(i*dF_dX_a+dF_dX)/(i+1);
  }

  if (flag==ON) {
    outputfile=efopen(outputfilename,"w");
    fprintf(outputfile,"%8.5lf ",dF_dX_a);
  }
  else {
    outputfile=efopen(outputfilename,"a");
    fprintf(outputfile,"%8.5lf \n",dF_dX_a);
  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename outputfilename parmfilename na1 na2 \n",progname);
}
