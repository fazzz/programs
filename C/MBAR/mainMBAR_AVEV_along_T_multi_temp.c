
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <getopt.h>

#include "EF.h"
#include "IO.h"

#include "MBAR.h"

#define ON 1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,t;
  int iflag=OFF,jflag=OFF;
  int interval=1;
  int d=0,num,numk;
  double sum=0.0;
  double f;

  double din;

  int n_sim;
  int *n;

  double T,Tmin,Tmax,dT,*T_sim;
  double k_B=1.98723e-3;
  double beta;

  double *fene;
  double ***enek,**U;
  double **od_data;

  double *AVEV;

  double criteria_BAR=1.0e-7;
  int MAXITE=1000;

  char *progname;
  char *enekfilename,*eneklistfilename;
  char *Tsimlistfilename;
  char *datalistfilename;
  char *outfilename;

  FILE *enekfile,*eneklistfile;
  FILE *Tsimlistfile;
  FILE *datalistfile;
  FILE *datafile;
  FILE *outfile;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  int opt_idx=1;

  struct option long_opt[] = {
    {"interval",1,NULL,'p'},
    {"Tmin",1,NULL,'m'},
    {"Tmax",1,NULL,'a'},
    {"dT",1,NULL,'d'},
    {"i",0,NULL,'i'},
    {"j",0,NULL,'j'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  progname=argv[0];
  while((c=getopt_long(argc,argv,"hip:m:a:d:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'p':
      interval=atoi(optarg);
      break;
    case 'm':
      Tmin=atof(optarg);
      break;
    case 'a':
      Tmax=atof(optarg);
      break;
    case 'd':
      dT=atof(optarg);
      break;
    case 'i':
      iflag=ON;
      break;
    case 'j':
      jflag=ON;
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  n_sim = atoi(*argv);
  eneklistfilename = *++argv;
  Tsimlistfilename = *++argv;
  datalistfilename = *++argv;
  outfilename = *++argv;

  fene=(double *)gcemalloc(sizeof(double)*n_sim);

  n=(int *)gcemalloc(sizeof(int)*n_sim);
  enek=(double ***)gcemalloc(sizeof(double **)*n_sim);
  U=(double **)gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<n_sim;++i) U[i]=(double *)gcemalloc(sizeof(double )*1);
  for (i=0;i<n_sim;++i) enek[i]=(double **)gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<n_sim;++i) enek[i]=(double **)gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<n_sim;++i) for (j=0;j<n_sim;++j) enek[i][j]=(double *)gcemalloc(sizeof(double ));
  od_data=(double **)gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<n_sim;++i) od_data[i]=(double **)gcemalloc(sizeof(double ));
  T_sim=(double *)gcemalloc(sizeof(double)*n_sim);

  Tsimlistfile=efopen(Tsimlistfilename,"r");
  for (i=0;i<n_sim;++i) fscanf(Tsimlistfile,"%lf",&T_sim[i]);
  fclose(Tsimlistfile);

  eneklistfile=efopen(eneklistfilename,"r");
  for (i=0;i<n_sim;++i) {
    getline(&line,&len,eneklistfile);
    l=strlen(line);
    line[l-1]='\0';
    enekfile=efopen(line,"r");
    getline(&line,&len,enekfile);
    k=0;
    num=0;
    d = 1;
    while ( d != -1  )  {
      d=fscanf(enekfile,"%d",&l);
      d=fscanf(enekfile,"%lf",&f);
      if (k%interval == 0) {
	U[i]=(double *)gcerealloc(U[i],sizeof(double)*(num+1));
	U[i][num]=f;
	++num;
      }
      ++k;
    } 
    fclose(enekfile);
    n[i]=num-1;
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
      if (jflag==ON)
	d=fscanf(datafile,"%d",&l);
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

  AVEV=(double *)gcemalloc(sizeof(double)*n_sim);
  for (i=0;i<n_sim;++i) {
    T=T_sim[i];
    for (j=0;j<n_sim;++j) {
      beta=1.0/k_B*(1.0/T_sim[j]-1.0/T);
      for (k=0;k<n_sim;++k) {
	num=0;
	for (l=0;l<n[k];++l) {
	  f=U[k][l];
	  enek[j][k]=(double *)gcerealloc(enek[j][k],sizeof(double)*(num+1));
	  enek[j][k][num]=exp(-1.0*f*beta);
	  ++num;
	} 
      }
    }

    MBAR_ite_high_speed_2(fene,enek,n_sim,n,criteria_BAR,MAXITE);
    AVEV[i]=MBAR_AVEV_multi_temp(fene,enek,n_sim,n,od_data,i);
  }

  outfile=efopen(outfilename,"w");
  fprintf(outfile,"T AVEV  \n");  
  for (i=0;i<n_sim;++i) {
    T=T_sim[i];
    fprintf(outfile,"%10.4lf %10.4lf\n",T,AVEV[i]);
  }
  fclose(outfile);

  return 0;
}

void USAGE(char *progname) {
  printf("-p -- interval\n");
  printf("-h -- help\n");
  printf("USAGE: %s n_sim width fenefilename enelistfilename datalistfilename histfilename histfilename pmffilename errorfile perrorfilename \n", progname);
}

