#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "FFL.h"

#include "TACCM.h"

#include "PTL.h"
#include "EF.h"
#include "NC.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numatom,numres,numstep,numZ;

  double *Z;

  int addflag=OFF;

  int nobetaflag=OFF;

  int AMBERMODEflag=OFF;

  double T=300,T_sim=300;
  double k_B=1.98723e-3;
  double beta=1.0;

  struct potential e;
  struct force f;
  double p_t=0.0;

  double KZ=10.0;

  int *indexTACCM,**pairsZ;
  double *theta,*f2;
  char *TACCMfilename,*Ztrjfilename;
  FILE *TACCMfile,*Ztrjfile;

  double x[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MD;
  struct my_netcdf_out_id_AMBER nc_id;

  double pi;

  char *inputfilename,*Thetafilename,*reffilename,*outputfilename,*parmfilename;
  FILE *Thetafile,*reffile,*outputfile,*parmfile;

  char *progname;

  int opt_idx=1;

  pi=acos(-1.0);

  struct option long_opt[] = {
    {"Amber",0,NULL,'A'},
    {"tempobj",1,NULL,'t'},
    {"tempsim",1,NULL,'s'},
    {"a",0,NULL,'a'},
    {"nobeta",0,NULL,'n'},
    {"KZ",1,NULL,'K'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hnaAt:K:s:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'A':
      AMBERMODEflag=ON;
      break;
    case 't':
      T=atof(optarg);
      break;
    case 's':
      T_sim=atof(optarg);
      break;
    case 'a':
      addflag=ON;
      break;
    case 'n':
      nobetaflag=ON;
      break;
    case 'K':
      KZ=atof(optarg);
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

  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }
  inputfilename     = *argv;
  Thetafilename     = *++argv;
  Ztrjfilename      = *++argv;
  parmfilename      = *++argv;
  TACCMfilename     = *++argv;
  outputfilename    = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  numres=AP.NRES;

  if (AMBERMODEflag==ON) numstep=mync_get_present_step_AMBER(inputfilename,&nc_id);
  else numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MD);

  TACCMfile=efopen(TACCMfilename,"r");
  fscanf(TACCMfile,"%d",&numZ);
  pairsZ=(int **)gcemalloc(sizeof(int *)*numZ);
  for (i=0;i<numZ;++i) pairsZ[i]=(int *)gcemalloc(sizeof(int)*5);
  for (i=0;i<numZ;++i) {
    for (j=0;j<4;++j) 
      fscanf(TACCMfile,"%d",&pairsZ[i][j]);
    fscanf(TACCMfile,"%d",&pairsZ[i][j]);
  }
  fclose(TACCMfile);
  theta=(double *)gcemalloc(sizeof(double)*numZ);
  f2=(double *)gcemalloc(sizeof(double)*numZ);
  Z=(double *)gcemalloc(sizeof(double)*numZ);

  for (i=0;i<numZ;++i) Z[i]=theta[i];

  ffL_set_calcffandforce(&e,&f);

  Thetafile=efopen(Thetafilename,"r");
  Ztrjfile=efopen(Ztrjfilename,"r");
  if ( nobetaflag==OFF )
    beta=1.0/k_B*(1.0/T_sim-1.0/T);
  if ( addflag==OFF ) outputfile=efopen(outputfilename,"w");
  else if ( addflag==ON )  outputfile=efopen(outputfilename,"a");
  for (i=0;i<numstep;++i) {
    for (j=0;j<numZ;++j) fscanf(Thetafile,"%lf",&theta[j]);
    for (j=0;j<numZ;++j) fscanf(Ztrjfile,"%lf",&Z[j]);

    p_t=TACCM_calc_eff_FF_Z(Z,numZ,theta,KZ,f2,pi);
    fprintf(outputfile,"%d %e \n",i,p_t*beta);
  }
  fclose(outputfile);
  fclose(Thetafile);
  fclose(Ztrjfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] Thetafilename Ztrjfilename parmfilename TACCMfilename outputfilename\n",progname);
}


