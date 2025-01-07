
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

  double *fel;

  double pi,minx,maxx,miny,maxy,dx,dy;

  double minxD, maxxD, minyD, maxyD;

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
  
  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }
  minxD = atof(*argv);
  maxxD = atof(*++argv);
  if (minxD>=maxxD) {
    printf("error\n");
    exit(1);
  }
  minyD = atof(*++argv);
  maxyD = atof(*++argv);
  if (minyD>=maxyD) {
    printf("error\n");
    exit(1);
  }
  inputfilename  = *++argv;
  outputfilename = *++argv;

  inputfile=efopen(inputfilename,"r");

  fel=(double *)gcemalloc(sizeof(double)*1);
  i=0;
  Nx=1;
  Ny=0;
  flag=OFF;
  while ((c=fscanf(inputfile,"%lf ",&fx)) != EOF) {
    fscanf(inputfile,"%lf ",&fy);

    if ((fx>= minxD && fx< maxxD) && (fy>= minyD && fy< maxyD)) {
      if (i==0) minx=fx;

      if ((fx>= minxD && fx< maxxD) && (fy>= minyD && fy< maxyD)) {
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
      }
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
    
    if ((fx>= minxD && fx< maxxD) && (fy>= minyD && fy< maxyD)) {
      if (flagnum==ON)
	fel[i]=fp;
      else 
	fel[i]=-1.0;
      ++i;
      fel=(double *)gcerealloc(fel,sizeof(double)*i);
    }

    if ((fx>= minxD && fx< maxxD) && (fy>= minyD && fy< maxyD)) {
      fx_old=fx;
      fy_old=fy;
    }
    maxx=fx_old;
  }
  dx=(maxx-minx)/(Nx-1);
  dy=(maxy-miny)/(Ny-1);

  fclose(inputfile);

  outputfile=efopen(outputfilename,"w");
  n=0;
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {

      x[0]=minx+dx*i;
      x[1]=miny+dy*j;

      if (fel[n]!=-1.0)
	fprintf(outputfile,"%12.8lf %12.8lf %12.8lf\n",x[0],x[1],fel[n]);
      else 
	fprintf(outputfile,"%12.8lf %12.8lf ******\n",x[0],x[1]);
      ++n;
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
