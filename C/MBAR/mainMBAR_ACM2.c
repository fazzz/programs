
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
  int interval=1;
  int d=0,num;
  double f;
  int n_sim;
  int *n;

  double *fene;
  double ***enek;
  double *covm;

  char *progname;
  char *enekfilename,*eneklistfilename,*fenefilename;
  char *outfilename;

  FILE *enekfile,*eneklistfile,*fenefile;
  FILE *logfile,*outfile;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"hp:"))!=-1) {
    switch(c) {
    case 'p':
      interval=atoi(optarg);
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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  n_sim = atoi(*argv);
  fenefilename = *++argv;
  eneklistfilename = *++argv;
  outfilename = *++argv;

  fene=(double *)gcemalloc(sizeof(double)*n_sim);
  fenefile=efopen(fenefilename,"r");
  getline(&line,&len,fenefile);
  for (i=0;i<n_sim;++i)
    fscanf(fenefile,"%d %lf",&d,&fene[i]);
  fclose(fenefile);

  n=(int *)gcemalloc(sizeof(int)*n_sim);
  enek=(double ***)gcemalloc(sizeof(double **)*n_sim);
  for (i=0;i<n_sim;++i) {
    enek[i]=(double **)gcemalloc(sizeof(double *)*n_sim);
  }
  for (i=0;i<n_sim;++i) {
    for (j=0;j<n_sim;++j) {
      enek[i][j]=(double *)gcemalloc(sizeof(double));
    }
  }

  eneklistfile=efopen(eneklistfilename,"r");
  for (i=0;i<n_sim;++i) {
    for (j=0;j<n_sim;++j) {
      getline(&line,&len,eneklistfile);
      line[strlen(line)-1]='\0';
      enekfile=efopen(line,"r");
      getline(&line,&len,enekfile);
      k=0;
      num=0;
      d = 1;
      while ( d != -1  )  {
	d=fscanf(enekfile,"%lf",&f);
	d=fscanf(enekfile,"%lf",&f);
	if (k%interval == 0) {
	  enek[i][j]=(double *)gcerealloc(enek[i][j],sizeof(double)*(num+1));
	  enek[i][j][num]=f;
	  ++num;
	}
	++k;
      } 
      fclose(enekfile);
      n[i]=num-1;
    }
  }
  fclose(eneklistfile);

  covm=MBAR_ACM2(fene,enek,n_sim,n);

  logfile=efopen("log_cov.txt","w");
  for (i=0;i<n_sim;++i)
    fprintf(logfile,"%e ",covm[i*n_sim+i]-2.0*covm[i*n_sim+0]+covm[0]);
  fprintf(logfile,"\n");
  fclose(logfile);

  outfile=efopen(outfilename,"w");
  fprintf(outfile,"# cov \n");
  for (i=0;i<n_sim;++i) {
    for (j=0;j<n_sim;++j)
      fprintf(outfile,"%e ",covm[i*n_sim+j]);
    fprintf(outfile,"\n");
  }
  fclose(outfile);
  
  return 0;
}

void USAGE(char *progname) {
  printf("-p -- interval\n");
  printf("-h -- help\n");
  printf("USAGE: %s n_sim width fenefilename enelistfilename outfilename \n", progname);
}

