
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "mymath.h"
#include "netcdf_mine.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int HEADERflag=OFF,HEADEROUTflag=OFF,NETCDFflag=OFF;
  int numrow=1;
  int interval=1;
  int nts,nbts=10,nmts=1,numpoints;
  int d;
  double f;
  double **timeseries,**timeseriesba;

  char *progname;
  char *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"HLnhr:b:m:"))!=-1) {
    switch(c) {
    case 'H':
      HEADERflag=ON;
      break;
    case 'L':
      HEADEROUTflag=ON;
      break;
    case 'n':
      NETCDFflag=ON;
      break;
    case 'r':
      numrow=atoi(optarg);
      break;
    case 'b':
      nbts=atoi(optarg);
      break;
    case 'm':
      nmts=atoi(optarg);
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

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  outputfilename = *++argv;

  timeseries=(double *)gcemalloc(sizeof(double)*numrow*1);

  inputfile=efopen(inputfilename,"r");
  nts=0;
  if (HEADERflag==ON)  {
    getline(&line,&len,inputfile);
    while ( d != -1  )  {
      for (j=0;j<numrow;++j) {
	d=fscanf(inputfile,"%d",&l);
	d=fscanf(inputfile,"%lf",&f);
	if (k%interval == 0) {
	  timeseries[j]=(double *)gcerealloc(timeseries[j],sizeof(double)*(nts+1));
	  timeseries[j][nts]=f;
	  ++nts;
	}
      }
      ++k;
    }
  }
  else {
    while ( d != -1  )  {
      for (j=0;j<numrow;++j) {
	d=fscanf(inputfile,"%lf",&f);
	if (k%interval == 0) {
	  timeseries[j]=(double *)gcerealloc(timeseries[j],sizeof(double)*(nts+1));
	  timeseries[j][nts]=f;
	  ++nts;
	}
      }
      ++k;
    }
  }
  fclose(inputfile);

  if (nts < nbts || nbts < nmts) {
    printf("error\n");
    exit(1);
  }

  numpoints=(int)nts/nmts;
  timeseriesba=(double **)gcemalloc(sizeof(double *)*numrow);
  for (i=0;i<numrow;++i) timeseriesba[i]=(double *)gcemalloc(sizeof(double)*1);
  numpoints=ts_bl_ave(timeseries,nts,numrow,nbts,nmts,timeseriesba);

  outputfile=efopen(outputfilename,"w");
  if (HEADEROUTflag==OFF) {
    for (i=0;i<numpoints;++i) {
      for (j=0;j<numrow;++j) fprintf(outputfile,"%12.8e ",timeseriesba[j][i]);
      fprintf(outputfile,"\n ");
    }
  }
  else {
    fprintf(outputfile,"# temp \n ");
    for (i=0;i<numpoints;++i) {
      fprintf(outputfile,"%d ",i+1);
      for (j=0;j<numrow;++j) fprintf(outputfile,"%12.8e ",timeseriesba[j][i]);
      fprintf(outputfile,"\n ");
    }
  }

    fclose(outputfile);

  return 0;
}

void USAGE(char *progname) {
  printf("USAGE: \n");
  printf("[-H] header option \n");
  printf("%s inputfilename(data) outputfilename(data_norm)\n",progname);
}
