
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EF.h"

#define ON 0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,n;

  int N;
  int Nx,Ny;
  int numx,numy;

  double minx,miny,maxx,maxy;
  double fx,fx_old,fy,fy_old;
  double dx,dy;

  int flag,flagnum,floatflag,signflag,expflag;
  int dim;

  double *pmf1,*error_pmf1,**pmf2,**error_pmf2,pmf_old;
  double **path,**path2,**error_path;

  double length,dist;

  double fp[2],ex,signex;

  int numpoint;

  char *inputfilename1,*inputfilename2,*outputfilename;
  FILE *inputfile1,*inputfile2,*outputfile;

  double pi;
  
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
      USAGE(progname);   exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  numpoint = atoi(*argv);
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  outputfilename = *++argv;

  inputfile1=efopen(inputfilename1,"r");
  pmf1=(double *)gcemalloc(sizeof(double)*1);
  error_pmf1=(double *)gcemalloc(sizeof(double)*1);
  i=0;
  Nx=1;
  Ny=0;
  flag=OFF;
  signflag=OFF;
  expflag=OFF;
  while ((c=fscanf(inputfile1,"%lf ",&fx)) != EOF) {
    if (i==0) minx=fx;
    
    fscanf(inputfile1,"%lf ",&fy);

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

    fp[0]=0.0;
    fp[1]=0.0;
    n=0;
    ex=0.0;
    floatflag=OFF;
    flagnum=OFF;
    while ((c=fgetc(inputfile1))!='\n') {
      if (c >= '0' && c <= '9') {
	flagnum=ON;
	if (expflag==ON)
	  ex=ex*10.0+c-'0';
	else if (floatflag==OFF && expflag==OFF)
	  fp[n]=fp[n]*10.0+c-'0';
	else if (floatflag==ON && expflag==OFF) {
	  dim+=1;
	  fp[n]+=(c-'0')/pow(10.0,dim);
	}
      }
      else if (c=='.') {
	floatflag=ON;
	dim=0;
      }
      else if (c=='e') {
	signflag=ON;
      }
      else if (c=='+') {
	if (signflag==ON) {
	  signex=1.0;
	  expflag=ON;
	}	  
      }
      else if (c=='-') {
	if (signflag==ON) {
	  signex=-1.0;
	  expflag=ON;
	}	  
      }
      else if (c==' ') {
	if (n==0) {
	  if (expflag==ON) {
	    ex=signex*ex;
	    fp[n]=fp[n]*pow(10.0,ex);
	    ex=0.0;
	    expflag=OFF;
	    signflag=OFF;
	  }

	  floatflag=OFF;
	  n=1;
	}
      }
    }     

    if (flagnum==ON) {
      pmf1[i]=fp[0];
      error_pmf1[i]=fp[1];
    }
    else {
      pmf1[i]=-1.0;
      error_pmf1[i]=-1.0;
    }

    fx_old=fx;
    fy_old=fy;
    ++i;
    pmf1=(double *)gcerealloc(pmf1,sizeof(double)*i);
    error_pmf1=(double *)gcerealloc(error_pmf1,sizeof(double)*i);
  }
  maxx=fx;

  dx=(maxx-minx)/(Nx-1);
  dy=(maxy-miny)/(Ny-1);

  fclose(inputfile1);

  pmf2=(double **)gcemalloc(sizeof(double *)*Nx);
  error_pmf2=(double **)gcemalloc(sizeof(double *)*Nx);
  for (i=0;i<Nx;++i) {
    pmf2[i]=(double *)gcemalloc(sizeof(double)*Ny);
    error_pmf2[i]=(double *)gcemalloc(sizeof(double)*Ny);
  }

  n=0;
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {
      pmf2[i][j]=pmf1[n];
      error_pmf2[i][j]=error_pmf1[n];
      ++n;
    }
  }

  path=(double **)gcemalloc(sizeof(double *)*numpoint);
  for (i=0;i<numpoint;++i) path[i]=(double *)gcemalloc(sizeof(double)*2);
  inputfile2=efopen(inputfilename2,"r");
  for (i=0;i<numpoint;++i) {
    fscanf(inputfile2,"%lf %lf",&path[i][0],&path[i][1]);
  }
  fclose(inputfile2);

  length=sqrt((path[numpoint-1][0]-path[0][0])*(path[numpoint-1][0]-path[0][0])+(path[numpoint-1][1]-path[0][1])*(path[numpoint-1][1]-path[0][1]));

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numpoint;++i) {
    numx=(int)((path[i][0]-minx)/dx);
    numy=(int)((path[i][1]-miny)/dy);

    dist=sqrt((path[i][0]-path[0][0])*(path[i][0]-path[0][0])+(path[i][1]-path[0][1])*(path[i][1]-path[0][1]));

    if ( pmf2[numx][numy] != pmf_old && pmf2[numx][numy] != -1.0 )
      fprintf(outputfile,"%8.3lf %8.3lf %8.3lf\n",dist/length,pmf2[numx][numy],error_pmf2[numx][numy]);

    pmf_old=pmf2[numx][numy];

  }
  fclose(outputfile);

  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename \n",progname);
}

