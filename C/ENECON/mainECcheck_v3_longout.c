#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ENECON_2.h"
#include "readAOUT.h"
#include "LEASQDV.h"

#include "EF.h"

#define ON 1
#define OFF 0

int usage(void);

int main(int argc, char *argv[]) {
  int i,j;
  int numatom,numstep=100,numini=0;
  int numrow=1,specrow=1;

  int moderd=ON,modead=ON,models=ON;
  int logflag=OFF,numsflag=OFF;

  double *ene;
  double rmsd,adiv,a,b;
  double f;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *progname;
  char *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;

  while((c=getopt(argc,argv,"haleLAn:i:o:s:"))!=-1) {
    switch(c) {
    case 'a':
      modead=OFF;
      break;
    case 'i':
      numini=atoi(optarg);
      break;
    case 'l':
      models=OFF;
      break;
    case 'L':
      logflag=ON;
      break;
    case 'o':
      numrow=atoi(optarg);
      break;
    case 's':
      specrow=atoi(optarg);
      break;
    case 'n':
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
  
  for (i=0;i<numini;++i) getline(&line,&len,inputfile);

  ene=(double *)gcemalloc(sizeof(double)*(numstep-numini+1));
  for (i=numini;i<numstep;++i) {
    for (j=0;j<numrow;++j) {
      fscanf(inputfile,"%lf",&f);
      if (j==specrow-1) {
	ene[i-numini]=f;
      }
    }
  }

  fclose(inputfile);

  if (moderd==ON) rmsd=ec_flu_2(ene,numstep-numini);
  if (modead==ON) adiv=ec_avd_2(ene,numstep-numstep);
  if (models==ON) least_sqrt_devi(ene,numstep-numini,&a,&b);

  if (logflag==ON) {
    if (rmsd!=0.0) rmsd=log10(rmsd);
    if (adiv!=0.0) adiv=log10(adiv);
    if (a!=0.0) a=log10(fabs(a));
  }


  outputfile=efopen(outputfilename,"w");
  if (logflag!=ON) {
    if (moderd==ON) fprintf(outputfile,"%24.20e\n",rmsd);
    if (modead==ON) fprintf(outputfile,"%24.20e\n",adiv);
    if (models==ON) fprintf(outputfile,"%24.20e %24.20e\n",a,b);
  }
  else {
    if (moderd==ON) fprintf(outputfile,"%24.20lf\n",rmsd);
    if (modead==ON) fprintf(outputfile,"%24.20lf\n",adiv);
    if (models==ON) fprintf(outputfile,"%24.20lf %24.20lf\n",a,b);
  }


  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-L ] log10  \n");
  printf("[-h ] help \n");
  printf("-l %s inputfilename outputfilename \n",progname);
}
