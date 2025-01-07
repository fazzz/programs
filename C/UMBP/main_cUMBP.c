#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "UMBP.h"

#include "PTL.h"
#include "EF.h"
#include "NC.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numatom,numstep;

  double *crd;
  double *pUMB,pUMB_t,**fUMB;

  int AMBERMODEflag=OFF;

  int *pairp,nump;
  double *fcp,*dihe_equp;

  double x[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id;

  double pi;

  char *inputfilename,*parmfilename,*umbfilename;
  char *trjfilename,*outputfilename;

  FILE *inputfile,*parmfile;
  FILE *outputfile;
  FILE *umbfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"Amber",0,NULL,'A'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hA",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'A':
      AMBERMODEflag=ON;
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  pi=acos(-1.0);

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  inputfilename     = *argv;
  parmfilename      = *++argv;
  umbfilename       = *++argv;
  outputfilename    = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  umbfile=efopen(umbfilename,"r");
  fscanf(umbfile,"%d",&nump);
  pairp=(int *)gcemalloc(sizeof(int)*nump*4);
  fcp=(double *)gcemalloc(sizeof(double)*nump);
  dihe_equp=(double *)gcemalloc(sizeof(double)*nump);
  pUMB=(int *)gcemalloc(sizeof(int)*nump);
  fUMB=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i)
    fUMB[i]=(double *)gcemalloc(sizeof(double)*3);
  for (i=0;i<nump;++i)
    fscanf(umbfile,"%d %d %d %d %lf %lf",&pairp[i*4+0],&pairp[i*4+1],&pairp[i*4+2],&pairp[i*4+3],&fcp[i],&dihe_equp[i]);
  fclose(umbfile);
  for (i=0;i<nump;++i) {
    for (j=0;j<4;++j) pairp[i*4+j]-=1;
    dihe_equp[i]=dihe_equp[i]*pi/180.0;
  }

  if (AMBERMODEflag==ON) numstep=mync_get_present_step_AMBER(inputfilename,&nc_id);
  else numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MCD);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    if (AMBERMODEflag==ON) mync_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id,crd_nc);
    else mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
    for (j=0;j<numatom;++j)for (k=0;k<3;++k)crd[j*3+k]=crd_nc[j][k];

    pUMB_t=UMB_calc_dihetype_ff(crd,numatom,pairp,nump,fcp,dihe_equp,pUMB,fUMB);

    fprintf(outputfile,"%e\n",pUMB_t);
  }
  fclose(outputfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parmfilename umbfilename outputfilename\n",progname);
}


