
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "Gaussian.h"
#include "RAND.h"
#include "EF.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,n;

  int n_Gaussian;
  int n_Gaussian_1/*=97*/;
  int n_Gaussian_2=5;
  int n_Gaussian_3=3;

  int nx,ny,a,b;

  double d;

  double x[2],minx=-8.0, maxx=-1.5, dx=0.2, miny=-9.0, maxy=8.0, dy=0.2,dx2,dy2;
  double **nyu, ***Sigma, *pi_k;
  double **nyu2, ***Sigma2, *pi_k2;
  double *prob,*prob2,fel,fel2;

  double pi,sum=0.0;

  char *outputfilename;
  FILE *outputfile;

  char *line;
  size_t len=0;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;
  
  pi=acos(-1.0);
  
  struct option long_opt[] = {
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"h",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);  exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < 1) {
    USAGE(progname);
    exit(1);
  }
  outputfilename = *argv;

  //  n_Gaussian=n_Gaussian_1+n_Gaussian_2;

  nx=(int)((maxx-minx)/dx/2);
  ny=(int)((maxy-miny)/dy/2);

  n_Gaussian_1=nx*ny;

  n_Gaussian=n_Gaussian_1+n_Gaussian_2;

  nyu=(double **)gcemalloc(sizeof(double *)*n_Gaussian);
  Sigma=(double ***)gcemalloc(sizeof(double **)*n_Gaussian);
  pi_k=(double *)gcemalloc(sizeof(double)*n_Gaussian);

  nyu2=(double **)gcemalloc(sizeof(double *)*n_Gaussian_3);
  Sigma2=(double ***)gcemalloc(sizeof(double **)*n_Gaussian_3);
  pi_k2=(double *)gcemalloc(sizeof(double)*n_Gaussian_3);

  for (i=0;i<n_Gaussian;++i) {
    nyu[i]=(double *)gcemalloc(sizeof(double)*2);

    Sigma[i]=(double **)gcemalloc(sizeof(double *)*2);
    for (j=0;j<2;++j) Sigma[i][j]=(double *)gcemalloc(sizeof(double )*2);
  }

  for (i=0;i<n_Gaussian_3;++i) {
    nyu2[i]=(double *)gcemalloc(sizeof(double)*2);

    Sigma2[i]=(double **)gcemalloc(sizeof(double *)*2);
    for (j=0;j<2;++j) Sigma2[i][j]=(double *)gcemalloc(sizeof(double )*2);
  }

  //  dx2=(maxx-minx)/n_Gaussian_2;
  //  dy2=(maxy-miny)/n_Gaussian_2;

  i=0;
  //  for (i=0;i<n_Gaussian_1;++i) {
  for (a=0;a<nx;++a) {
    for (b=0;b<ny;++b) {
      nyu[i]=(double *)gcemalloc(sizeof(double)*2);
      nyu[i][0]=minx+dx*a*2/*(i-n_Gaussian)*/;
      nyu[i][1]=miny+dy*b*2/*(i-n_Gaussian)*/;

      Sigma[i]=(double **)gcemalloc(sizeof(double *)*2);
      for (j=0;j<2;++j) Sigma[i][j]=(double *)gcemalloc(sizeof(double )*2);
      d=genrand_real2();
      Sigma[i][0][0]=0.001*d;
      Sigma[i][0][1]=0.00001*d;
      Sigma[i][1][0]=0.00001*d;
      Sigma[i][1][1]=0.001*d;

      pi_k[i]=/*0.000001*genrand_real2()*/0.0;

      /************************************************/
      /* Sigma[i][0][0]=genrand_real2()*0.1;	      */
      /* Sigma[i][0][1]=0.0/\*genrand_real2()*0.5*\/; */
      /* Sigma[i][1][0]=0.0/\*genrand_real2()*0.5*\/; */
      /* Sigma[i][1][1]=genrand_real2()*0.1;	      */
      /* 					      */
      /* pi_k[i]=0.0/\*genrand_real2()*0.1*\/;	      */
      /************************************************/
      //    sum+=pi_k[i];
      ++i;
    }
  }

  /****************************************/
  /* for (i=0;i<n_Gaussian_1;++i) {	  */
  /*   nyu[i][0]=minx+dx2*i;		  */
  /*   nyu[i][1]=miny+dy2*i;		  */
  /* 					  */
  /*   Sigma[n_Gaussian_1-1][0][0]=0.01;  */
  /*   Sigma[n_Gaussian_1-1][0][1]=0.001; */
  /*   Sigma[n_Gaussian_1-1][1][0]=0.004; */
  /*   Sigma[n_Gaussian_1-1][1][1]=0.02;  */
  /* 					  */
  /*   pi_k[i]=1.0;			  */
  /* }					  */
  /****************************************/

  nyu[n_Gaussian_1][0]=-5.0;
  nyu[n_Gaussian_1][1]=-5.0;

  nyu[n_Gaussian_1+1][0]=-6.0;
  nyu[n_Gaussian_1+1][1]=6.0;
  
  nyu[n_Gaussian_1+2][0]=-4.0;
  nyu[n_Gaussian_1+2][1]=-2.0;
  
  nyu[n_Gaussian_1+3][0]=-3.0;
  nyu[n_Gaussian_1+3][1]=4.0;
  
  nyu[n_Gaussian_1+4][0]=-4.5;
  nyu[n_Gaussian_1+4][1]=4.0;

  nyu2[0][0]=-5.5;
  nyu2[0][1]=1.0;

  nyu2[1][0]=-6.0;
  nyu2[1][1]=2.0;
  
  nyu2[2][0]=-5.0;
  nyu2[2][1]=-1.0;

  Sigma[n_Gaussian_1][0][0]=0.005;
  Sigma[n_Gaussian_1][0][1]=0.0;
  Sigma[n_Gaussian_1][1][0]=0.0;
  Sigma[n_Gaussian_1][1][1]=0.005;

  Sigma[n_Gaussian_1+1][0][0]=0.0001;
  Sigma[n_Gaussian_1+1][0][1]=0.0;
  Sigma[n_Gaussian_1+1][1][0]=0.0;
  Sigma[n_Gaussian_1+1][1][1]=0.001;
  
  Sigma[n_Gaussian_1+2][0][0]=0.003;
  Sigma[n_Gaussian_1+2][0][1]=0.0;
  Sigma[n_Gaussian_1+2][1][0]=0.0;
  Sigma[n_Gaussian_1+2][1][1]=0.003;
  
  Sigma[n_Gaussian_1+3][0][0]=0.0008;
  Sigma[n_Gaussian_1+3][0][1]=0.0;
  Sigma[n_Gaussian_1+3][1][0]=0.0;
  Sigma[n_Gaussian_1+3][1][1]=0.0008;
  
  Sigma[n_Gaussian_1+4][0][0]=0.0005;
  Sigma[n_Gaussian_1+4][0][1]=0.0;
  Sigma[n_Gaussian_1+4][1][0]=0.0;
  Sigma[n_Gaussian_1+4][1][1]=0.0005;

  Sigma2[0][0][0]=0.004;
  Sigma2[0][0][1]=0.0;
  Sigma2[0][1][0]=0.0;
  Sigma2[0][1][1]=0.004;
  
  Sigma2[1][0][0]=0.0004;
  Sigma2[1][0][1]=0.0;
  Sigma2[1][1][0]=0.0;
  Sigma2[1][1][1]=0.0004;
  
  Sigma2[2][0][0]=0.002;
  Sigma2[2][0][1]=0.0;
  Sigma2[2][1][0]=0.0;
  Sigma2[2][1][1]=0.002;

  pi_k[n_Gaussian_1]=1000.0;
  pi_k[n_Gaussian_1+1]=60.0;
  pi_k[n_Gaussian_1+2]=100.0;
  pi_k[n_Gaussian_1+3]=100.0;
  pi_k[n_Gaussian_1+4]=-100.0;

  pi_k2[0]=1.0;
  pi_k2[1]=1.0;
  pi_k2[2]=1.0;

  /********************************/
  /* sum=0.0;			  */
  /* for (i=0;i<n_Gaussian;++i) { */
  /*   sum+=pi_k[i];		  */
  /* }				  */
  /* for (i=0;i<n_Gaussian;++i) { */
  /*   pi_k[i]=pi_k[i]/sum;	  */
  /* }				  */
  /********************************/

  prob=Create_mixed_twoD_GaussianMap(minx, maxx, dx,
				     miny, maxy, dy,
				     nyu,  Sigma, 
				     pi_k, n_Gaussian, pi);

  prob2=Create_mixed_twoD_GaussianMap(minx, maxx, dx,
				      miny, maxy, dy,
				      nyu2,  Sigma2, 
				      pi_k2, n_Gaussian_3, pi);

  outputfile=efopen(outputfilename,"w");
  n=0;
  x[0]=minx;
  for (i=0;x[0]<maxx;++i) {
    x[0]=minx+dx*i;
    x[1]=miny;
    for (j=0;x[1]<maxy;++j){
      x[1]=miny+dy*j;
      fel=-1.0*log(prob[n]);
      if (fel>1000.0) fel=1000.0;
      if (fel<0.0) fel=0.0;
      fel2=-1.0*log(prob2[n]);
      if (fel2>1000.0) fel2=1000.0;
      if (fel2<0.0) fel2=0.0;
      fprintf(outputfile,"%8.3lf %8.3lf %8.3lf\n",x[0],x[1],fel-0.2*fel2+genrand_real2()*30.0);
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
