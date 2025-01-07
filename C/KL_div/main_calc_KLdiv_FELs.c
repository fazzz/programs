
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
  int i,j,k,l,n;
  double fx,fx_old,fy_old,fy,fp;

  int index_minx,index_maxx,index_miny,index_maxy;

  double pi,x_map[2],minx, maxx, dx=0.01, miny, maxy, dy=0.01,minpmf;

  double minxD, maxxD, minyD, maxyD;

  int N,K;
  int Nx,Ny;
  int Nx_old,Ny_old;

  int dim;
  int flag,floatflag,flagnum,periodflag=OFF,flagd=OFF,flagminmaxcheck[4];

  double **prob,*prob2,sum;
  double **prob_app,*prob2_app,sum_app;

  double KL;

  char *inputfilename_pmf,*inputfilename_apmf;
  FILE *inputfile_pmf,*inputfile_apmf;
  
  char *line;
  size_t len=0;
  
  int c,d;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;
  
  //  pi=acos(-1.0);
  
  struct option long_opt[] = {
    {"d",0,NULL,'d'},
    {"minx",1,NULL,'x'},
    {"maxx",1,NULL,'a'},
    {"miny",1,NULL,'y'},
    {"maxy",1,NULL,'m'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"hdx:a:y:m:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'd':
      flagd=ON;  break;
    case 'x':
      minxD=atof(optarg);  break;
    case 'a':
      maxxD=atof(optarg);  break;
    case 'y':
      minyD=atof(optarg);  break;
    case 'm':
      maxyD=atof(optarg);  break;
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
  inputfilename_pmf  = *argv;
  inputfilename_apmf = *++argv;

  inputfile_pmf=efopen(inputfilename_pmf,"r");

  prob2=(double *)gcemalloc(sizeof(double)*1);
  i=0;
  Nx=1;
  Ny=0;
  flag=OFF;
  while ((c=fscanf(inputfile_pmf,"%lf ",&fx)) != EOF) {
    fscanf(inputfile_pmf,"%lf ",&fy);

    if (flagd == OFF) {
      if (i==0) minx=fx;
    }
    else {
      if ((fx>= minxD && fx< maxxD) && (fy>= minyD && fy< maxyD)) {
	if (i==0) minx=fx;
      }
    }

    if (flagd == OFF) {
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
    else {
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
    while ((c=fgetc(inputfile_pmf))!='\n') {
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
	dim=0/*1*//*2012-10-17*/;
      }
      else {
	flagnum=OFF;
      }
    }     

    if (flagd==OFF) {
      if (flagnum==ON)
	//      if (fp <10)
	prob2[i]=exp(-1.0*fp);
      /**************************/
      /* else		        */
      /* 	prob2[i]=0.0;   */
      /**************************/
      else 
	prob2[i]=0.0;
      ++i;
      prob2=(double *)gcerealloc(prob2,sizeof(double)*i);
    }
    else {
      if ((fx>= minxD && fx< maxxD) && (fy>= minyD && fy< maxyD)) {
	if (flagnum==ON)
	  //      if (fp <10)
	  prob2[i]=exp(-1.0*fp);
	/**************************/
	/* else		        */
	/* 	prob2[i]=0.0;   */
	/**************************/
	else 
	  prob2[i]=0.0;
	++i;
	prob2=(double *)gcerealloc(prob2,sizeof(double)*i);
      }
    }

    if (flagd==OFF) {
      fx_old=fx;
      fy_old=fy;
    }
    else {
      if ((fx>= minxD && fx< maxxD) && (fy>= minyD && fy< maxyD)) {
	fx_old=fx;
	fy_old=fy;
      }
    }
  }
  if (flagd==OFF)
    maxx=fx;
  else
    maxx=fx_old;

  dx=(maxx-minx)/(Nx-1);
  dy=(maxy-miny)/(Ny-1);

  fclose(inputfile_pmf);

  inputfile_apmf=efopen(inputfilename_apmf,"r");

  prob2_app=(double *)gcemalloc(sizeof(double)*1);
  i=0;
  Nx=1;
  Ny=0;
  flag=OFF;
  while ((c=fscanf(inputfile_apmf,"%lf ",&fx)) != EOF) {
    fscanf(inputfile_apmf,"%lf ",&fy);

    if (flagd == OFF) {
      if (i==0) minx=fx;
    }
    else {
      if ((fx>= minxD && fx< maxxD) && (fy>= minyD && fy< maxyD)) {
	if (i==0) minx=fx;
      }
    }

    if (flagd == OFF) {
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
    else {
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
    while ((c=fgetc(inputfile_apmf))!='\n') {
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
	dim=0/*1*//*2012-10-17*/;
      }
      else {
	flagnum=OFF;
      }
    }     

    if (flagd==OFF) {
      if (flagnum==ON)
	//      if (fp <10)
	prob2_app[i]=exp(-1.0*fp);
      /**************************/
      /* else		        */
      /* 	prob2[i]=0.0;   */
      /**************************/
      else 
	prob2_app[i]=0.0;
      ++i;
      prob2_app=(double *)gcerealloc(prob2_app,sizeof(double)*i);
    }
    else {
      if ((fx>= minxD && fx< maxxD) && (fy>= minyD && fy< maxyD)) {
	if (flagnum==ON)
	  //      if (fp <10)
	  prob2_app[i]=exp(-1.0*fp);
	/**************************/
	/* else		        */
	/* 	prob2[i]=0.0;   */
	/**************************/
	else 
	  prob2_app[i]=0.0;
	++i;
	prob2_app=(double *)gcerealloc(prob2_app,sizeof(double)*i);
      }
    }

    if (flagd==OFF) {
      fx_old=fx;
      fy_old=fy;
    }
    else {
      if ((fx>= minxD && fx< maxxD) && (fy>= minyD && fy< maxyD)) {
	fx_old=fx;
	fy_old=fy;
      }
    }
  }
  if (flagd==OFF)
    maxx=fx;
  else
    maxx=fx_old;

  dx=(maxx-minx)/(Nx-1);
  dy=(maxy-miny)/(Ny-1);

  fclose(inputfile_apmf);

  prob=(double **)gcemalloc(sizeof(double *)*Nx);
  prob_app=(double **)gcemalloc(sizeof(double *)*Nx);
  for (i=0;i<Nx;++i) {
    prob[i]=(double *)gcemalloc(sizeof(double)*Ny);
    prob_app[i]=(double *)gcemalloc(sizeof(double)*Ny);
  }

  n=0;
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {
      prob[i][j]=prob2[n];
      prob_app[i][j]=prob2_app[n];
      ++n;
    }
  }

  sum=0.0;
  sum_app=0.0;
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {
      sum+=prob[i][j];
      sum_app+=prob_app[i][j];
    }
  }
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {
      prob[i][j]=prob[i][j]/sum;
      prob_app[i][j]=prob_app[i][j]/sum_app;
    }
  }

  KL=0.0;
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {
      KL+=prob[i][j]*log(prob[i][j]/prob_app[i][j]);
    }
  }

  KL=KL/(Nx*Ny);

  printf("%12.8lf\n",KL);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename_pmf inputfilename_apmf\n",progname);
}
