#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

#include "PTL.h"
#include "EF.h"
#include "NC.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numatom,numres,numstep;

  double *crd,*refcrd;

  int addflag=OFF;

  double ep=ep_natatt_PROTEINS2008;

  int nibnum=3,criteria=criteria_PROTEINS2008;

  double T=300;
  double k_B=1.98723e-3;
  double KBT;

  struct potential e;
  struct force f;
  struct potential_GOLMAA_PROTEINS2008 e_GOLM;

  double x[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id;

  char *inputfilename,*reffilename,*outputfilename,*parmfilename;
  FILE *inputfile,*reffile,*outputfile,*parmfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"temp",1,NULL,'t'},
    {"a",0,NULL,'a'},
    {"ep",1,NULL,'e'},
    {"cutoff",1,NULL,'c'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hat:e:c:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 't':
      T=atof(optarg);
      break;
    case 'a':
      addflag=ON;
      break;
    case 'e':
      ep=atof(optarg);
      break;
    case 'c':
      criteria=atof(optarg);
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
  numres=AP.NRES;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  reffile=efopen(reffilename,"r");
  getline(&line,&len,reffile);
  fscanf(reffile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(reffile,"%lf",&refcrd[i*3+j]);
  fclose(reffile);

  numstep=mync_get_present_step_MCD(inputfilename,&nc_id);

  ffL_set_calcffandforce(&e,&f);
  GOLMAA_PROTEINS2008_ff_set_calcff(&e_GOLM,refcrd,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);

  KBT=k_B*T;
  if ( addflag==OFF ) outputfile=efopen(outputfilename,"w");
  else if ( addflag==ON )  outputfile=efopen(outputfilename,"a");
  for (i=0;i<numstep;++i) {
    mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id,crd_nc);
    for (j=0;j<numatom;++j)
      for (k=0;k<3;++k) crd[j*3+k]=crd_nc[j][k];

    GOLMAA_PROTEINS2008_ff_calcff(crd,numatom,&e_GOLM);
    fprintf(outputfile,"%d %e \n",i,e_GOLM.p_t/KBT);
  }
  fclose(outputfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfilename parmfilename outputfilename\n",progname);
}


