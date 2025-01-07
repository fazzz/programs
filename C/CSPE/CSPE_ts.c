
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>

#include "fftw3.h"

#include "IO.h"
#include "EF.h"
#include "SPE.h"
#include "netcdf_mine.h"

#define ON 1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,n;

  int vflag=OFF,nflag=OFF;

  int numstep,num;
  int num_of_split=1,numini=1,numfin,numinirow,numfinrow,numrow;

  double c=2.999792e-2;
  double kb=1.98723e-3*4.18407*100.0;

  double dt,pi;
  double **spe,**dat,**speave;

  char *inputfilename,*outputfilename,*progname;
  FILE *inputfile,*outputfile;

  int ncid,indexofsim_dimid,fene_varid;
  size_t start[1],count[1];

  char *line,*dummy;
  size_t len=0;

  int c2;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c2=getopt(argc,argv,"hnvi:f:j:g:d:s:r:"))!=-1) {
    switch(c2) {
    case 'n':
      nflag=ON;
      break;
    case 'v':
      vflag=ON;
      break;
    case 'i':
      numini=atof(optarg)-1;
      break;
    case 'f':
      numfin=atoi(optarg)-1;
      break;
    case 'j':
      numinirow=atoi(optarg)-1;
      break;
    case 'g':
      numfinrow=atoi(optarg)-1;
      break;
    case 'd':
      dt=atof(optarg);
      break;
    case 's':
      num_of_split=atoi(optarg);
      break;
    case 'r':
      numrow=atoi(optarg);
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

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  outputfilename = *++argv;

  inputfile=efopen(inputfilename,"r");
  if (nflag==OFF) for (i=0;i<numini;++i)  getline(&line,&len,inputfile);

  numstep=(int)((numfin-numini+1)/num_of_split);
  num=numfinrow-numinirow+1;
  speave=(double **)gcemalloc(sizeof(double *)*num);
  for (i=0;i<num;++i) speave[i]=(double *)gcemalloc(sizeof(double)*numstep);
  for (i=0;i<num_of_split;++i) {
    spe=(double **)gcemalloc(sizeof(double *)*num);
    for (j=0;j<num;++j) spe[j]=(double *)gcemalloc(sizeof(double)*numstep);
    dat=(double **)gcemalloc(sizeof(double *)*num);
    for (j=0;j<num;++j) dat[j]=(double *)gcemalloc(sizeof(double)*numstep);
    if (nflag==ON) {
      ;
    }
    else io_scanmatrix_ifv(inputfile,numstep,numrow,numinirow,numfinrow,dat);
    for (j=0;j<num;++j) FT_ts(numstep,dat[j],spe[j]);
    for (j=0;j<numstep;++j) for (k=0;k<num;++k)	speave[k][j]=(i*speave[k][j]+spe[k][j])/(i+1);
  }
  fclose(inputfile);
    
  outputfile=efopen(outputfilename,"w");
  pi=acos(-1.0);
  for (i=0;i<numstep;++i) {
    fprintf(outputfile,"%e ",(double)i/numstep/dt/c);
    for (j=0;j<num;++j)
      if (vflag==ON) fprintf(outputfile,"%e ",speave[j][i]);
      else fprintf(outputfile,"%e ",(2.0*pi*i/numstep/dt)*(2.0*pi*i/numstep/dt)*speave[j][i]);
    fprintf(outputfile,"\n ");
  }
  fclose(outputfile);

  return 0;
}

void USAGE(char *progname) {
  printf("-n -- netcdf output  \n");
  printf("-v -- derivative mode  \n");
  printf("-i -- numini\n");
  printf("-f -- numfinp\n");
  printf("-j -- numinirow\n");
  printf("-g -- numfinrow\n");
  printf("-d -- deltat\n");
  printf("-s -- num of split s\n");
  printf("-r -- numrow\n");
  printf("-h -- help\n");
  printf("USAGE:%s [-n] [-v] [-i] [-j] [-g] [-d] [-s] [-r] [-h] inputfilename outputfile\n", progname);
}
