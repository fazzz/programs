#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "GOLMAA.h"
#include "GOLMAA_set.h"
#include "GOLMAA_check.h"
#include "GOLMAA_dbasin.h"
#include "PTL.h"
#include "EF.h"

#include "netcdf_mine.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int initialstep,numstep,numatom,numatomp;
  double ratio;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crd,*refcrd1,*refcrd2;
  double **crdA,**crdB;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id_MD;

  double delta=0.5,deltaV=0.5;

  struct potential_GOLMAA_dbasin e;
  double R_C_D=1.0;

  char *inputfilename,*reffilename1,*reffilename2,*parmfilename,*outputfilename;
  char *progname;
  FILE *inputfile,*parmfile,*reffile1,*reffile2,*outputfile;

  while((c=getopt(argc,argv,"haAKi:d:V:"))!=-1) {
    switch(c) {
    case 'i':
      initialstep=atoi(optarg);
      break;
    case 'd':
      delta=atof(optarg);
      break;
    case 'V':
      deltaV=atof(optarg);
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

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  reffilename1 = *++argv;
  reffilename2 = *++argv;
  parmfilename = *++argv;
  outputfilename = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd1=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd2=(double *)gcemalloc(sizeof(double)*numatom*3);

  reffile1=efopen(reffilename1,"r");
  io_scanconf_Amber_ini(reffile1,numatom,refcrd1);
  fclose(reffile1);

  reffile2=efopen(reffilename2,"r");
  io_scanconf_Amber_ini(reffile2,numatom,refcrd2);
  fclose(reffile2);

  numstep=mync_get_present_step_AMBER(inputfilename,&nc_id_MD);

  GOLMAAff_dbasin_set_calcff(&e,refcrd1,refcrd2,numatom,R_C_D);
  outputfile=efopen(outputfilename,"w");

  for (i=0;i<numstep;++i) {
    mync_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id_MD,crd_nc);
    for (j=0;j<numatom;++j)
      for (k=0;k<3;++k) 
	crd[j*3+k]=crd_nc[j][k];

    ratio=GOLMAAff_dbasin_calcff_ratio(crd,numatom,&e,delta,deltaV);

    fprintf(outputfile,"%lf\n",ratio);

  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-b] CA \n");
  printf("[-c] HV \n");
  printf("[-h] help \n");
  printf("%s [-b] [-c] [-h] inputfilename reffilename1 reffilename2 parmfilename outputfilename \n",progname);
}
