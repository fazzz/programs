
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>

#include "EF.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,n;
  double f;

  int numstep,num;
  int num_of_split=1,numini=1,numfin,numinirow,numfinrow,numrow;

  double c=2.999792e-2;
  double kb=1.98723e-3*4.18407*100.0;

  double dt,sum=0.0;

  int numfiles;

  char **inputfilenames,*outputfilename,*progname;
  FILE **inputfiles,*outputfile;

  char *line,*dummy;
  size_t len=0;

  int c2;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c2=getopt(argc,argv,"hoi:f:j:g:d:s:r:"))!=-1) {
    switch(c2) {
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

  if (argc < 1) {
    USAGE(progname);
    exit(1);
  }
  numfiles       = atoi(*argv);
  inputfilenames = (char **)gcemalloc(sizeof(char *)*numfiles);
  inputfiles     = (FILE **)gcemalloc(sizeof(FILE *)*numfiles);
  for (i=0;i<numfiles;++i) inputfilenames[i]  = *++argv;
  outputfilename = *++argv;

  numstep=(int)((numfin-numini+1)/num_of_split);
  num=numfinrow-numinirow/*+1*/;
  if (num<1) {
    printf("error: numfinrow must be larger or equal to numinirow\n");
    exit(1);
  }

  for (i=0;i<numfiles;++i)
    inputfiles[i]=efopen(inputfilenames[i],"r");
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {

    for (j=0;j<numfiles;++j) {
      fscanf(inputfiles[j],"%lf",&f);
      fscanf(inputfiles[j],"%lf",&f);
    }

    fprintf(outputfile,"%e %e ",(double)i/numstep/dt/c,(double)numstep*dt/i);

    for (j=0;j<num;++j) {
      sum=0.0;
      for (k=0;k<numfiles;++k) {
	fscanf(inputfiles[k],"%lf",&f);
	sum+=f;
      }
      sum=sum/numfiles;

      fprintf(outputfile,"%e ",sum);
    }

    fprintf(outputfile,"\n");
  }
  for (i=0;i<numfiles;++i) fclose(inputfiles[i]);
  
  return 0;
}

void USAGE(char *progname) {
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

