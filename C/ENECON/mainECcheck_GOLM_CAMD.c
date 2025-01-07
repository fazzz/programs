#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "ENECON.h"
#include "readAOUT.h"
#include "LEASQDV.h"

#include "EF.h"

#include "netcdf_mine.h"

#define ON 1
#define OFF 0

int usage(void);

int main(int argc, char *argv[]) {
  int i,j;
  int numatom,numstep=1;

  int moderd=ON,modead=ON,models=ON;
  int logflag=OFF;

  double *PE,*KE,*KEv,*PEv,*ET,*T,*E_NC;
  double rmsd,adiv,a,b;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id;

  char *progname;
  char *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;

  while((c=getopt(argc,argv,"hraleLAn:N:"))!=-1) {
    switch(c) {
    case 'r':
      moderd=OFF;
      break;
    case 'a':
      modead=OFF;
      break;
    case 'l':
      models=OFF;
      break;
    case 'L':
      logflag=ON;
      break;
    case 'n':
      numstep=atoi(optarg);
      break;
    case 'N':
      inputfilename=optarg;
      numstep=mync_get_present_step_MCD(inputfilename,&nc_id);
      break;
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  outputfilename = *++argv;

  PE=(double *)gcemalloc(sizeof(double)*numstep);
  KE=(double *)gcemalloc(sizeof(double)*numstep);
  KEv=(double *)gcemalloc(sizeof(double)*numstep);
  PEv=(double *)gcemalloc(sizeof(double)*numstep);
  PE=(double *)gcemalloc(sizeof(double)*numstep);
  ET=(double *)gcemalloc(sizeof(double)*numstep);
  T=(double *)gcemalloc(sizeof(double)*numstep);
  E_NC=(double *)gcemalloc(sizeof(double)*numstep);

  inputfile=efopen(inputfilename,"r");
  for (i=0;i<numstep;++i) {
    fscanf(inputfile,"%d",&j);
    fscanf(inputfile,"%lf",&PE[i]);
    fscanf(inputfile,"%lf",&KE[i]);
    fscanf(inputfile,"%lf",&KEv[i]);
    fscanf(inputfile,"%lf",&PEv[i]);
    fscanf(inputfile,"%lf",&ET[i]);
    fscanf(inputfile,"%lf",&T[i]);
    fscanf(inputfile,"%lf",&E_NC[i]);
  }
  fclose(inputfile);

  if (moderd==ON) rmsd=ec_flu(ET,numstep);
  if (modead==ON) adiv=ec_avd(ET,numstep);
  if (models==ON) least_sqrt_devi(ET,numstep,&a,&b);

  if (logflag==ON) {
    if (rmsd!=0.0) rmsd=log10(rmsd);
    if (adiv!=0.0) adiv=log10(adiv);
    if (a!=0.0) a=log10(fabs(a));
  }

  outputfile=efopen(outputfilename,"w");
  if (moderd==ON) fprintf(outputfile,"%14.10lf\n",rmsd);
  if (modead==ON) fprintf(outputfile,"%14.10lf\n",adiv);
  if (models==ON) fprintf(outputfile,"%14.10lf %14.10lf\n",a,b);
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-L ] log10  \n");
  printf("[-n ] numstep  \n");
  printf("[-N ] neccdffilename  \n");
  printf("[-h ] help \n");
  printf("-l %s inputfilename outputfilename \n",progname);
}
