
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "IO.h"

#include "WHAM.h"
#include "MBAR.h"

#define ON 1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,t;
  int flag=OFF;
  int interval=1;
  int d=0,num;
  double sum=0.0;
  double f;

  double din;

  int n_sim;
  int *n;

  double *fene;
  double ***enek;

  int frame;
  double max,min;
  double width;
  double **od_data,*hist,*pmf;
  double pmf_min;

  int n_total=0;
  double *W,*WT;
  double *hist_error,*pmf_error;

  double criteria_BAR=1.0e-7;
  int MAXITE=1000;

  char *progname;
  char *enekfilename,*eneklistfilename,*fenefilename;
  char *datalistfilename,*datafilename;
  char *histfilename,*pmffilename;
  char *errorfilename,*perrorfilename;

  FILE *enekfile,*eneklistfile,*fenefile;
  FILE *datalistfile,*datafile;
  FILE *histfile,*pmffile;
  FILE *errorfile,*perrorfile;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"hmp:"))!=-1) {
    switch(c) {
    case 'm':
      flag=ON;
      break;
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

  if (argc < 7) {
    USAGE(progname);
    exit(1);
  }
  n_sim = atoi(*argv);
  width = atof(*++argv);
  fenefilename = *++argv;
  eneklistfilename = *++argv;
  datalistfilename = *++argv;
  histfilename = *++argv;
  pmffilename = *++argv;
  errorfilename = *++argv;
  perrorfilename = *++argv;

  fene=(double *)gcemalloc(sizeof(double)*n_sim);
  if (flag==OFF) {
    fenefile=efopen(fenefilename,"r");
    getline(&line,&len,fenefile);
    for (i=0;i<n_sim;++i)
      fscanf(fenefile,"%d %lf",&d,&fene[i]);
    fclose(fenefile);
  }

  n=(int *)gcemalloc(sizeof(int)*n_sim);
  enek=(double ***)gcemalloc(sizeof(double **)*n_sim);
  for (i=0;i<n_sim;++i) {
    enek[i]=(double **)gcemalloc(sizeof(double *)*n_sim);
  }
  for (i=0;i<n_sim;++i) {
    for (j=0;j<n_sim;++j) {
      enek[i][j]=(double *)gcemalloc(sizeof(double ));
    }
  }
  od_data=(double **)gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<n_sim;++i) {
    od_data[i]=(double **)gcemalloc(sizeof(double ));
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
	  enek[j][i]=(double *)gcerealloc(enek[j][i],sizeof(double)*(num+1));
	  enek[j][i][num]=f;
	  ++num;
	}
	++k;
      } 
      fclose(enekfile);
      n[i]=num-1;
    }
  }
  fclose(eneklistfile);

  datalistfile=efopen(datalistfilename,"r");
  for (i=0;i<n_sim;++i) {
    getline(&line,&len,datalistfile);
    line[strlen(line)-1]='\0';
    datafile=efopen(line,"r");
    k=0;
    num=0;
    d = 1;
    while ( d != -1  )  {
      d=fscanf(datafile,"%lf",&f);
      if (k%interval == 0) {
	od_data[i]=(double *)gcerealloc(od_data[i],sizeof(double)*(num+1));
	od_data[i][num]=f;
	++num;
      }
      ++k;
    } 
    fclose(datafile);
  }
  fclose(datalistfile);

  if (flag==ON) WHAM_fast_BFGS(fene,enek,n_sim,n);
  for (i=0;i<n_sim;++i) n_total+=n[i];
    
  hist=MBAR_AVE_oned(fene,enek,n_sim,n,od_data,width,&max,&min,&frame);

  sum=0.0;
  for (i=0;i<frame;++i) {
    sum+=hist[i];
  }
  for (i=0;i<frame;++i) {
    hist_error[i]=hist_error[i]/sum;
  }

  pmf=(double *)gcemalloc(sizeof(double)*frame);
  pmf_error=(double *)gcemalloc(sizeof(double)*frame);
  pmf_min=0.0;
  for (i=0;i<=frame;++i) {
    if ( pmf_min == 0.0 && hist[i] > 0.0 )
	pmf_min=hist[i];
    if ( pmf_min > pmf[i] && hist[i] > 0.0 )
	pmf_min=hist[i];
  }


  for (i=0;i<=frame;++i) {
    if (hist[i]!=0.0) {
      pmf[i]=-log(hist[i])+log(pmf_min);
    }
    if (hist_error[i]!=0.0) {
      pmf_error[i]=-log(hist_error[i]);
    }
  }

  if (flag==ON) {
    fenefile=efopen(fenefilename,"w");
    fprintf(fenefile,"# fene\n");
    for (i=0;i<n_sim;++i)
      fprintf(fenefile,"%d %12.10lf\n",i,fene[i]);
    fclose(fenefile);
  }

  histfile=efopen(histfilename,"w");
  fprintf(histfile,"# hist \n");
  for (i=0;i<frame;++i)
    fprintf(histfile,"%10.4lf %10.4lf\n",min+width*i,hist[i]);
  fclose(histfile);

  pmffile=efopen(pmffilename,"w");
  fprintf(pmffile,"# pmf \n");
  for (i=0;i<frame;++i)
    fprintf(pmffile,"%10.4lf %10.4lf\n",min+width*i,pmf[i]);
  fclose(pmffile);

  return 0;
}

void USAGE(char *progname) {
  printf("-p -- interval\n");
  printf("-h -- help\n");
  printf("USAGE: %s n_sim width fenefilename enelistfilename datalistfilename histfilename histfilename pmffilename errorfile perrorfilename \n", progname);
}

