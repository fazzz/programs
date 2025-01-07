#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PTLb.h"
#include "FFLc.h"

#include "EF.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,m;
  int numatomv,numatomw,numatomw_w;

  int d;
  
  char *line;
  size_t len=0;

  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crdv,*crdw,*crdw_w;

  struct potential ev,ew,ew_w;
  struct force fv,fw,fw_w;
  struct AmberParmL apv,apw,apw_w;

  int c;

  char *crdvfilename,*crdwfilename,*crdw_wfilename,*parmfilenamev,*parmfilenamew,*parmfilenamew_w,*progname;
  FILE *crdvfile,*crdwfile,*crdw_wfile,*parmfilev,*parmfilew,*parmfilew_w;

  int opt_idx=1;

  struct option long_opt[] = {
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"h",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }
  crdvfilename    = *argv;
  crdwfilename    = *++argv;
  crdw_wfilename  = *++argv;
  parmfilenamev   = *++argv;
  parmfilenamew   = *++argv;
  parmfilenamew_w = *++argv;

  parmfilev=efopen(parmfilenamev,"r");
  readParmtopLb(parmfilev,&apv);
  fclose(parmfilev);
  numatomv=apv.NATOM;

  parmfilew=efopen(parmfilenamew,"r");
  readParmtopLb(parmfilew,&apw);
  fclose(parmfilew);
  numatomw=apw.NATOM;

  parmfilew_w=efopen(parmfilenamew_w,"r");
  readParmtopLb(parmfilew_w,&apw_w);
  fclose(parmfilew_w);
  numatomw_w=apw_w.NATOM;

  crdv=(double *)gcemalloc(sizeof(double)*numatomv*3);
  crdw=(double *)gcemalloc(sizeof(double)*numatomw*3);
  crdw_w=(double *)gcemalloc(sizeof(double)*numatomw_w*3);

  crdvfile=efopen(crdvfilename,"r");
  getline(&line,&len,crdvfile);
  fscanf(crdvfile,"%d",&d);
  for (i=0;i<numatomv;++i) for (j=0;j<3;++j) fscanf(crdvfile,"%lf",&crdv[i*3+j]);
  fclose(crdvfile);

  crdwfile=efopen(crdwfilename,"r");
  getline(&line,&len,crdwfile);
  fscanf(crdwfile,"%d",&d);
  for (i=0;i<numatomw;++i) for (j=0;j<3;++j) fscanf(crdwfile,"%lf",&crdw[i*3+j]);
  fclose(crdwfile);

  crdw_wfile=efopen(crdw_wfilename,"r");
  getline(&line,&len,crdw_wfile);
  fscanf(crdw_wfile,"%d",&d);
  for (i=0;i<numatomw_w;++i) for (j=0;j<3;++j) fscanf(crdw_wfile,"%lf",&crdw_w[i*3+j]);
  fclose(crdw_wfile);

  ffLc_set_calcffandforce(&ev,&fv,apv);
  ffLc_set_calcffandforce(&ew,&fw,apw);
  ffLc_set_calcffandforce(&ew_w,&fw_w,apw_w);

  ffLc_calcffandforce(crdv,numatomv,&ev,&fv,apv);
  ffLc_calcffandforce(crdw,numatomw,&ew,&fw,apw);
  ffLc_calcffandforce(crdw_w,numatomw_w,&ew_w,&fw_w,apw_w);

  printf("%10.6e \n",-1.0*(ew.p_t-ev.p_t-ew_w.p_t));
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [--Amber] [-h] crdvfilename crdwfilename crdw_wfilename parmtopvfilename parmtopwfilename parmtopw_wfilename\n",progname);
}
