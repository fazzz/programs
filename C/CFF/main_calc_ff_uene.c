#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "FFL.h"

#include "PTL.h"
#include "EF.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numatom,numres,numstep;

  double *crd;

  int addflag=OFF;
  int nobetaflag=OFF;

  int AMBERMODEflag=OFF;
  
  double T=300,T_sim=300;
  double k_B=1.98723e-3;
  double KBT;
  double beta=1.0;

  struct potential e;
  struct force f;

  double x[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MD;
  struct my_netcdf_out_id_AMBER nc_id;

  char *inputfilename,*outputfilename,*parmfilename;
  FILE *inputfile,*outputfile,*parmfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"Amber",0,NULL,'A'},
    {"temp",1,NULL,'t'},
    {"nobeta",0,NULL,'b'},
    {"a",0,NULL,'a'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hAabt:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'A':
      AMBERMODEflag=ON;
      break;
    case 't':
      T=atof(optarg);
      break;
    case 'a':
      addflag=ON;
      break;
    case 'b':
      nobetaflag=ON;
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
  inputfilename     = *argv;
  parmfilename      = *++argv;
  outputfilename    = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  numres=AP.NRES;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  if (AMBERMODEflag==ON) numstep=mync_get_present_step_AMBER(inputfilename,&nc_id);
  else numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MD);

  ffL_set_calcffandforce(&e,&f);
  ffL_calcffandforce(crd,numatom,&e,&f);

  if ( nobetaflag==OFF )
    beta=1.0/k_B*(1.0/T_sim-1.0/T);
  if ( addflag==OFF ) outputfile=efopen(outputfilename,"w");
  else if ( addflag==ON )  outputfile=efopen(outputfilename,"a");
  for (i=0;i<numstep;++i) {
    if (AMBERMODEflag==ON) mync_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id,crd_nc);
    else mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MD,crd_nc);
    for (j=0;j<numatom;++j)
      for (k=0;k<3;++k) crd[j*3+k]=crd_nc[j][k];

    ffL_calcffandforce(crd,numatom,&e,&f);
    fprintf(outputfile,"%d %e \n",i,e.p_t*beta);
  }
  fclose(outputfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfilename parmfilename outputfilename\n",progname);
}


