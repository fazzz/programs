#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "RAND.h"
#include "BOXMULL.h"

#include "EF.h"

#define ON  1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int uniformflag=OFF;
  int numini=100;
  double phi,psi,omega;
  double aphi=0.0,apsi=0.0,aomega=180.0;
  double vphi=30.0,vpsi=30.0,vomega=30.0;
  double maxphi=1.0,minphi=0.0;
  double maxpsi=1.0,minpsi=0.0;
  double maxomega=1.0,minomega=0.0;

  char *progname;
  char *angfilenamebase,angfilename[100];
  FILE *angfile;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  progname=argv[0];
  while((c=getopt(argc,argv,"hun:s:j:o:t:i:p:"))!=-1) {
    switch(c) {
    case 'u':
      uniformflag=ON;
      break;
    case 'n':
      numini=atoi(optarg);
      break;
    case 's':
      apsi=atof(optarg);
      maxpsi=atof(optarg);
      break;
    case 'j':
      aphi=atof(optarg);
      maxphi=atof(optarg);
      break;
    case 'o':
      aomega=atof(optarg);
      maxomega=atof(optarg);
      break;
    case 't':
      vpsi=atof(optarg);
      minpsi=atof(optarg);
      break;
    case 'i':
      vphi=atof(optarg);
      minphi=atof(optarg);
      break;
    case 'p':
      vomega=atof(optarg);
      minomega=atof(optarg);
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  argc-=optind;
  argv+=optind;

  if (argc < 1) {
    USAGE(progname);
    exit(1);
  }
  angfilenamebase  = *argv;

  for (i=0;i<numini;++i) {
    if (uniformflag==ON) {
      psi=(maxpsi-minpsi)*genrand_real2()+minpsi;
      phi=(maxphi-minphi)*genrand_real2()+minphi;
      omega=(maxomega-minomega)*genrand_real2()+minomega;
    }
    else {
      psi=Box_Muller(i,apsi,vpsi);
      phi=Box_Muller(i,aphi,vphi);
      omega=Box_Muller(i,aomega,vomega);
    }

    sprintf(angfilename,"%s_%d",angfilenamebase,i+1);
    angfile=efopen(angfilename,"w");
    for (j=0;j<10;++j)
      fprintf(angfile,"   0.000");
    fprintf(angfile,"\n");
    fprintf(angfile,"%8.3lf ",phi);
    fprintf(angfile,"%8.3lf ",psi);
    fprintf(angfile,"%8.3lf ",omega);
    for (j=0;j<7;++j)
      fprintf(angfile,"   0.000");
    fprintf(angfile,"\n");
    for (j=0;j<10;++j)
      fprintf(angfile,"   0.000");
    fprintf(angfile,"\n");

    fclose(angfile);
  }

  return 0;
}

void USAGE(char *progname) {
  printf("-n -- numini\n");
  printf("-u -- uniform distribution\n");
  printf("-s -- ave psi ( if u-mode max psi )\n");
  printf("-j -- ave phi ( if u-mode max phi )\n");
  printf("-o -- ave omega ( if u-mode max omega )\n");
  printf("-t -- variance psi ( if u-mode min psi )\n");
  printf("-i -- variance phi ( if u-mode min phi )\n");
  printf("-p -- variance omega ( if u-mode min omega )\n");
  printf("-h -- help\n");
  printf("USAGE: %s angfilename\n", progname);
}

