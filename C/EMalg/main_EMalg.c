
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EMalg.h"
#include "Gaussian.h"
#include "EF.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,n;
  double f;

  int N,K;

  double *prob;

  double **x_n,**nyu_k,***Sigma_k, *pi_k,**gamma_nk;
  double lnrou,lnrou_new,dlnrou,threhold=0.000001;

  double pi,x_map[2],minx, maxx, dx=0.01, miny, maxy, dy=0.01;

  char *inputfilename,*outputfilename,*logfilename;
  FILE *inputfile,*outputfile,*logfile;
  
  char *line;
  size_t len=0;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;
  
  pi=acos(-1.0);
  
  struct option long_opt[] = {
    {"dx",1,NULL,'x'},
    {"dy",1,NULL,'y'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"hx:y:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'x':
      dx=atof(optarg);  
      break;
    case 'y':
      dy=atof(optarg);  
      break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);   exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  K = atoi(*argv);
  inputfilename = *++argv;
  outputfilename = *++argv;

  inputfile=efopen(inputfilename,"r");
  i=0;
  x_n=(double **)gcemalloc(sizeof(double *)*1);
  x_n[0]=(double *)gcemalloc(sizeof(double)*2);
  while ((c=fscanf(inputfile,"%lf",&f))!=EOF) {
    x_n[i][0]=f;
    if (c=fscanf(inputfile,"%lf",&f)!=EOF)
      x_n[i][1]=f;
    else {
      exit(1);
      printf("inputfile error ! data not match");
    }
    x_n=(double **)gcerealloc(x_n,sizeof(double *)*(i+1));
    x_n[i+1]=(double *)gcemalloc(sizeof(double )*2);
    ++i;
  }
  fclose(inputfile);
  N=i;

  nyu_k=(double **)gcemalloc(sizeof(double *)*K);
  for (i=0;i<K;++i) nyu_k[i]=(double *)gcemalloc(sizeof(double)*2);

  Sigma_k=(double ***)gcemalloc(sizeof(double **)*K);
  for (i=0;i<K;++i) {
    Sigma_k[i]=(double **)gcemalloc(sizeof(double *)*2);
    for (j=0;j<2;++j) Sigma_k[i][j]=(double *)gcemalloc(sizeof(double )*2);
  }

  pi_k=(double *)gcemalloc(sizeof(double )*K);

  gamma_nk=(double **)gcemalloc(sizeof(double *)*N);
  for (i=0;i<N;++i) gamma_nk[i]=(double *)gcemalloc(sizeof(double)*K);

  for (i=0;i<K;++i) pi_k[i]=1.0;
  for (i=0;i<K;++i) {
    nyu_k[i][0]=x_n[i][0];
    nyu_k[i][1]=x_n[i][1];
  }
  for (i=0;i<K;++i) {
    Sigma_k[i][0][0]=1.0;
    Sigma_k[i][0][1]=0.0;
    Sigma_k[i][1][0]=0.0;
    Sigma_k[i][1][1]=1.0;
  }

  lnrou=EM_lnrou(N,K,x_n,nyu_k,Sigma_k,pi_k,pi);
  i=0;
  dlnrou=0.0;
  printf("# ofiteration = %3d ln(rou) = %8.3lf\n",i,lnrou);
  for (j=0;j<K;++j) {
    printf("nyu (%3d) = ( ",j);
    for (k=0;k<2;++k) {
      printf("%8.3lf ",nyu_k[j][k]);
    }
    printf(")\n");
  }

  while (dlnrou > threhold || i==0) {

    E_step(N,K,x_n,nyu_k,Sigma_k,pi_k,gamma_nk,pi);

    M_step(N,K,x_n,nyu_k,Sigma_k,pi_k,gamma_nk);
    
    lnrou_new=EM_lnrou(N,K,x_n,nyu_k,Sigma_k,pi_k,pi);

    dlnrou=fabs(lnrou_new-lnrou);
    lnrou=lnrou_new;

    ++i;
    printf("# ofiteration = %3d ln(rou) = %8.3lf\n",i,lnrou);
    for (j=0;j<K;++j) {
      printf("nyu (%3d) = ( ",j);
      for (k=0;k<2;++k) {
	printf("%8.3lf ",nyu_k[j][k]);
      }
      printf(")\n");
    }
  }

  minx=x_n[0][0];
  maxx=x_n[0][0];
  miny=x_n[0][1];
  maxy=x_n[0][1];
  for (i=0;i<N;++i) {
    if (minx>x_n[i][0]) minx=x_n[i][0];
    if (maxx<x_n[i][0]) maxx=x_n[i][0];
    if (miny>x_n[i][1]) miny=x_n[i][1];
    if (maxy<x_n[i][1]) maxy=x_n[i][1];
  }

  prob=Create_mixed_twoD_GaussianMap(minx, maxx, dx,
				     miny, maxy, dy,
				     nyu_k,  Sigma_k, pi_k, K,  pi);

  outputfile=efopen(outputfilename,"w");
  n=0;
  x_map[0]=minx;
  for (i=0;x_map[0]<maxx;++i) {
    x_map[0]=minx+dx*i;
    x_map[1]=miny;
    for (j=0;x_map[1]<maxy;++j){
      x_map[1]=miny+dy*j;
      fprintf(outputfile,"%8.3lf %8.3lf %10.6lf\n",x_map[0],x_map[1],prob[n]);
      ++n;
    }
  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename \n",progname);
}


