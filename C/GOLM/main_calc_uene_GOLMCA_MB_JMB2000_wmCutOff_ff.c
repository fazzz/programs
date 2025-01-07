#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "GOLM_Clementi_set_wmCutOff.h"
#include "GOLM_Clementi.h"
#include "GOLM_Clementi_MB_wmCutOff.h"

#include "PTL.h"
#include "EF.h"
#include "RAND.h"
#include "BOXMULL.h"
#include "MD.h"
#include "MD_NHC_MP1996.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,dummy;
  int numatom,numCAatom,numres,numstep;
  double pi;

  int addflag=OFF,nobetaflag=OFF;

  double ep=ep_natatt_Clementi;
  double de=1.0,d=1.0,d2;

  int nc=1;
  double T,T_sim,beta;
  double k_B=1.98723e-3;

  double cutoff=6.5;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MD;

  double *crd,*refcrd1,*refcrd2,*refcrdAA1,*refcrdAA2;

  struct potential e_AA;
  struct potential_GOLM_Clementi_MB e;

  double x[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*reffilename1,*reffilename2,*parmfilename;
  char *outputfilename;

  FILE *inputfile,*reffile1,*reffile2,*parmfile;
  FILE *outputfile,*rstfile,*rstvelfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"tempobj",1,NULL,'t'},
    {"tempsim",1,NULL,'s'},
    {"ep",1,NULL,'e'},
    {"de",1,NULL,'d'},
    {"dd",1,NULL,'2'},
    {"cutoff",1,NULL,'c'},
    {"a",0,NULL,'a'},
    {"nobeta",0,NULL,'n'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hant:s:e:d:2:c:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 't':
      T=atof(optarg);
      break;
    case 's':
      T_sim=atof(optarg);
      break;
    case 'e':
      ep=atof(optarg);
      break;
    case 'd':
      d=atof(optarg);
      break;
    case '2':
      de=atof(optarg);
      break;
    case 'c':
      cutoff=atof(optarg);
      break;
    case 'a':
      addflag=ON;
      break;
    case 'n':
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

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  inputfilename     = *argv;
  reffilename1      = *++argv;
  reffilename2      = *++argv;
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
  numres=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double)*numCAatom*3);
  refcrd1=(double *)gcemalloc(sizeof(double)*numCAatom*3);
  refcrd2=(double *)gcemalloc(sizeof(double)*numCAatom*3);
  refcrdAA1=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrdAA2=(double *)gcemalloc(sizeof(double)*numatom*3);

  reffile1=efopen(reffilename1,"r");
  getline(&line,&len,reffile1);
  fscanf(reffile1,"%d",&dummy);
  j=0;
  for (i=0;i<numatom;++i) {
    for (k=0;k<3;++k) fscanf(reffile1,"%lf",&refcrdAA1[i*3+k]);
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      for (k=0;k<3;++k) refcrd1[j*3+k]=refcrdAA1[i*3+k];
      ++j;
    }
  }
  fclose(reffile1);

  reffile2=efopen(reffilename2,"r");
  getline(&line,&len,reffile2);
  fscanf(reffile2,"%d",&dummy);
  j=0;
  for (i=0;i<numatom;++i) {
    for (k=0;k<3;++k) fscanf(reffile2,"%lf",&refcrdAA2[i*3+k]);
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      for (k=0;k<3;++k) refcrd2[j*3+k]=refcrdAA2[i*3+k];
      ++j;
    }
  }
  fclose(reffile2);

  numstep=myncL_get_present_step_MCD(inputfilename,&nc_id_MD);

  GOLM_Clementi_MB_ff_set_calcff_wmCutOff(&e,refcrd1,refcrd2,refcrdAA1,refcrdAA2,numCAatom,numatom,ep,cutoff);

  d2=d*d;
  d2=d2*k_B;
  de=de*k_B;

  GOLM_Clementi_MB_ff_calcff(crd,numCAatom,de,d2,&e);

  if ( nobetaflag==OFF )
    beta=1.0/k_B*(1.0/T_sim-1.0/T);
  if ( addflag==OFF ) outputfile=efopen(outputfilename,"w");
  else if ( addflag==ON )  outputfile=efopen(outputfilename,"a");

  for (i=0;i<numstep;++i) {
    myncL_open_inq_get_sh_MCD(inputfilename,numCAatom,i,1,i+1,&nc_id_MD,crd_nc);

    for (j=0;j<numCAatom;++j) for (k=0;k<3;++k) crd[j*3+k]=crd_nc[j][k];

    GOLM_Clementi_MB_ff_calcff(crd,numCAatom,de,d2,&e);
    fprintf(outputfile,"%d %e \n",i,e.p_MB*beta);
  }
  fclose(outputfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename(trajectory) refcrd1 refcrd2 parmfilename outputfilename\n",progname);
}
