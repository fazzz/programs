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
  int numatom,numstep,numini=0;

  int moderd=ON,modead=ON,models=ON;
  int NVE=OFF;
  int Amberflag=OFF;
  int logflag=OFF,numsflag=OFF;

  double *ene;
  double rmsd,adiv,a,b;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *progname;
  char *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;

  while((c=getopt(argc,argv,"hraleLAn:i:"))!=-1) {
    switch(c) {
    case 'r':
      moderd=OFF;
      break;
    case 'a':
      modead=OFF;
      break;
    case 'i':
      numini=atoi(optarg);
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

  inputfile=efopen(inputfilename,"r");

  if ( /*numstep*/numsflag==OFF && Amberflag==OFF && NVE==OFF)
    ene=readOUTwab(inputfile,"toal_vertial_energy      =",&numstep,numini);
  else if ( /*numstep*/numsflag==ON && Amberflag==OFF && NVE==OFF)
    ene=readOUTwabwnum(inputfile,"toal_vertial_energy      =",numstep,numini);
  if ( /*numstep*/numsflag==OFF && Amberflag==OFF && NVE==ON)
    ene=readOUTwab(inputfile,"toal_energy      =",&numstep,numini);
  else if ( /*numstep*/numsflag==ON && Amberflag==OFF && NVE==ON)
    ene=readOUTwabwnum(inputfile,"toal_energy      =",numstep,numini);                               
  else if ( /*numstep*/numsflag==OFF && Amberflag==ON )
    ene=readOUTwab(inputfile,"Etot   =  ",&numstep,numini);
  else if ( /*numstep*/numsflag==ON && Amberflag==ON )
    ene=readOUTwabwnum(inputfile,"Etot   =  ",numstep,numini);
  fclose(inputfile);

  if (moderd==ON) rmsd=ec_flu(ene,numstep);
  if (modead==ON) adiv=ec_avd(ene,numstep);
  if (models==ON) least_sqrt_devi(ene,numstep,&a,&b);

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
