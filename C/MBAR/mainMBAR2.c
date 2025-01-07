
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "IO.h"

#include "MBAR.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int d=0,num;
  int num_sim;
  int *numstep;
  double f;
  double ***ene;
  double **fene;

  char *progname;
  char *enefilename,*enelistfilename,*fenefilename;
  FILE *enefile,*enelistfile,*fenefile;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"h"))!=-1) {
    switch(c) {
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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  num_sim = atoi(*argv);
  enelistfilename  = *++argv;
  fenefilename = *++argv;

  numstep=(int *)gcemalloc(sizeof(int)*num_sim);
  ene=(double ***)gcemalloc(sizeof(double **)*num_sim);
  for (i=0;i<num_sim;++i) {
    ene[i]=(double **)gcemalloc(sizeof(double *)*num_sim);
  }
  for (i=0;i<num_sim;++i) {
    for (j=0;j<num_sim;++j) {
      ene[i][j]=(double *)gcemalloc(sizeof(double ));
    }
  }

  enelistfile=efopen(enelistfilename,"r");
  for (i=0;i<num_sim;++i) {
    for (j=0;j<num_sim;++j) {
      getline(&enefilename,&len,enelistfile);
      enefilename[strlen(enefilename)-1]='\0';
      enefile=efopen(enefilename,"r");
      getline(&line,&len,enefile);
      k=0;
      num=0;
      d = 1;
      while ( d != -1  )  {
	d=fscanf(enefile,"%lf",&f);
	d=fscanf(enefile,"%lf",&f);
	ene[i][j]=(double *)gcerealloc(ene[i][j],sizeof(double)*(k+1));
	ene[i][j][k]=f;
	++num;
	++k;
      } 
      fclose(enefile);
      numstep[i]=num-1;
    }
  }
  fclose(enelistfile);

  fene=(double **)gcemalloc(sizeof(double *)*num_sim);
  for (i=0;i<num_sim;++i)
    fene[i]=(double *)gcemalloc(sizeof(double)*numstep[i]);
  
  MBAR_ite2(fene,ene,num_sim,numstep);

  fenefile=efopen(fenefilename,"w");
  fprintf(fenefile,"# fene \n");
  for (i=0;i<num_sim;++i)
    fprintf(fenefile,"%d %10.4lf\n",i,fene[i]);
  fclose(fenefile);
  
  return 0;
}

void USAGE(char *progname) {
  printf("-h -- help\n");
  printf("USAGE: %s num_sim parmtopname enelistfilename fenefilename \n", progname);
}

