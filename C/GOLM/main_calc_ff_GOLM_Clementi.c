#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "GOLM_Clementi_set.h"
#include "GOLM_Clementi.h"

#include "PTL.h"
#include "EF.h"
#include "NC.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numatom,numCAatom,numstep;

  int addflag=OFF;

  double dt=0.001;
  double *crd,*refcrd,*refcrdAA;
  struct potential_GOLM_Clementi e;

  double T=300;
  double k_B=1.98723e-3;
  double KBT;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id;

  char *inputfilename,*parmfilename,*reffilename,*outputfilename;
  FILE *parmfile,*reffile,*outputfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"temp",1,NULL,'t'},
    {"a",0,NULL,'a'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hat:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 't':
      T=atof(optarg);
      break;
    case 'a':
      addflag=ON;
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  inputfilename     = *argv;
  reffilename       = *++argv;
  parmfilename      = *++argv;
  outputfilename    = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  j=0;
  for (i=0;i<numatom;++i) {
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      ++j;
    }
  }
  numCAatom=j;

  crd=(double *)gcemalloc(sizeof(double)*numCAatom*3);
  refcrd=(double *)gcemalloc(sizeof(double)*numCAatom*3);
  refcrdAA=(double *)gcemalloc(sizeof(double)*numatom*3);

  reffile=efopen(reffilename,"r");
  getline(&line,&len,reffile);
  fscanf(reffile,"%d",&d);
  j=0;
  for (i=0;i<numatom;++i) {
    for (k=0;k<3;++k) fscanf(reffile,"%lf",&refcrdAA[i*3+k]);
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      for (k=0;k<3;++k) refcrd[j*3+k]=refcrdAA[i*3+k];
      ++j;
    }
  }
  fclose(reffile);

  numstep=mync_get_present_step_MCD(inputfilename,&nc_id);

  GOLM_Clementi_ff_set_calcff(&e,refcrd,refcrdAA,numCAatom,numatom);

  KBT=k_B*T;
  if ( addflag==OFF ) outputfile=efopen(outputfilename,"w");
  else if ( addflag==ON )  outputfile=efopen(outputfilename,"a");
  for (i=0;i<numstep;++i) {
    mync_open_inq_get_sh_MCD(inputfilename,numCAatom,i,1,i+1,&nc_id,crd_nc);
    for (j=0;j<numCAatom;++j)
      for (k=0;k<3;++k) crd[j*3+k]=crd_nc[j][k];

    GOLM_Clementi_ff_calcff(crd,numCAatom,&e);
    fprintf(outputfile,"%d %e \n",i,e.p_t/KBT);
  }
  fclose(outputfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename reffilename parmfilename outputfilename\n",progname);
}
