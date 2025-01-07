#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ENECON.h"
#include "readAOUT.h"
#include "LEASQDV.h"

#include "EF.h"

#define ON 1
#define OFF 0

int usage(void);

int main(int argc, char *argv[]) {
  int i,j;
  int numatom,numstep=1;

  int moderd=ON,modead=ON,models=ON;
  int NVE=OFF;
  int Amberflag=OFF;
  int logflag=OFF,numsflag=OFF;

  double *PE,*KE,*KEv,*PEv,*PET,*s,*T,*Q_NC,dummy;
  double rmsd,adiv,a,b;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *progname;
  char *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;

  while((c=getopt(argc,argv,"hraleLAn:"))!=-1) {
    switch(c) {
    case 'r':
      moderd=OFF;
      break;
    case 'a':
      modead=OFF;
      break;
    case 'e':
      NVE=ON;
      break;
    case 'l':
      models=OFF;
      break;
    case 'L':
      logflag=ON;
      break;
    case 'A':
      Amberflag=ON;
      break;
    case 'n':
      numsflag=ON;
      numstep=atoi(optarg);
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
  PET=(double *)gcemalloc(sizeof(double)*numstep);
  s=(double *)gcemalloc(sizeof(double)*numstep);
  T=(double *)gcemalloc(sizeof(double)*numstep);
  Q_NC=(double *)gcemalloc(sizeof(double)*numstep);

  inputfile=efopen(inputfilename,"r");
  for (i=0;i<numstep;++i) {
    fscanf(inputfile,"%d",&j);
    fscanf(inputfile,"%lf",&PE[i]);
    fscanf(inputfile,"%lf",&KE[i]);
    fscanf(inputfile,"%lf",&KEv[i]);
    fscanf(inputfile,"%lf",&PEv[i]);
    fscanf(inputfile,"%lf",&PET[i]);
    fscanf(inputfile,"%lf",&s[i]);
    fscanf(inputfile,"%lf",&T[i]);
    fscanf(inputfile,"%lf",&Q_NC[i]);
    fscanf(inputfile,"%lf",&dummy);
  }
  fclose(inputfile);

  if (moderd==ON) rmsd=ec_flu(PET,numstep);
  if (modead==ON) adiv=ec_avd(PET,numstep);
  if (models==ON) least_sqrt_devi(PET,numstep,&a,&b);

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
  printf("[-h ] help \n");
  printf("-l %s inputfilename outputfilename \n",progname);
}
