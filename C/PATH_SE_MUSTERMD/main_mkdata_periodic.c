
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <netcdf.h>
#include <getopt.h>

#include "EF.h"

#define ON 0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,n,m;
  double fx,fx_old,fy_old,fy,fp;

  int N;
  int Nx,Ny;
  int Nx_old,Ny_old;

  int dim;
  int flag,floatflag,flagnum;

  double x[2];

  double **fel,*fel2;

  double pi,minx,maxx,miny,maxy,dx,dy;

  char *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;
  
  char *line;
  size_t len=0;
  
  int c,d;
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
  
  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  outputfilename = *++argv;

  inputfile=efopen(inputfilename,"r");

  fel2=(double *)gcemalloc(sizeof(double)*1);
  i=0;
  Nx=1;
  Ny=0;
  flag=OFF;
  while ((c=fscanf(inputfile,"%lf ",&fx)) != EOF) {
    fscanf(inputfile,"%lf ",&fy);

    if (i==0) minx=fx;

    if (fx != fx_old && i!=0) {
      if (flag==OFF)
	maxy=fy_old;
      ++Nx;
      flag=ON;
    }
    
    if (flag==OFF) {
      if (i==0)
	miny=fy;
      ++Ny;
    }

    fp=0.0;
    floatflag=OFF;
    while ((c=fgetc(inputfile))!='\n') {
      if (c >= '0' && c <= '9') {
	flagnum=ON;
	if (floatflag==OFF)
	  fp=fp*10.0+c-'0';
	else  {
	  dim+=1;
	  fp+=(c-'0')/pow(10.0,dim);
	}
      }
      else if (c=='.') {
	floatflag=ON;
	dim=0;
      }
      else {
	flagnum=OFF;
      }
    }     

    if (flagnum==ON)
      fel2[i]=fp;
    else 
      fel2[i]=0.0;
    ++i;
    fel2=(double *)gcerealloc(fel2,sizeof(double)*i);

    fx_old=fx;
    fy_old=fy;

  }
  maxx=fx;

  dx=(maxx-minx)/(Nx-1);
  dy=(maxy-miny)/(Ny-1);

  fclose(inputfile);

  fel=(double **)gcemalloc(sizeof(double *)*Nx*3);
  for (i=0;i<Nx*3;++i) {
    fel[i]=(double *)gcemalloc(sizeof(double)*Ny*3);
  }

  n=0;
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {

      fel[i][j]=fel2[n];

      fel[i+Nx][j]=fel2[n];
      fel[i+2*Nx][j]=fel2[n];

      fel[i][j+Ny]=fel2[n];
      fel[i][j+2*Ny]=fel2[n];

      fel[i+Nx][j+Ny]=fel2[n];
      fel[i+Nx][j+2*Ny]=fel2[n];

      fel[i+2*Nx][j+Ny]=fel2[n];
      fel[i+2*Nx][j+2*Ny]=fel2[n];

      ++n;
    }
  }

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<Nx*3;++i) {
    for (j=0;j<Ny*3;++j) {
      n=i-Nx;
      m=j-Ny;

      x[0]=minx+dx*n;
      x[1]=miny+dy*m;

      fprintf(outputfile,"%12.8lf %12.8lf %12.8lf\n",x[0],x[1],fel[i][j]);
    }
  }

  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename outputfilename \n",progname);
}
