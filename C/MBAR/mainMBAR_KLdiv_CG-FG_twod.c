
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "netcdf_mine.h"
#include "EF.h"
#include "IO.h"

#include "MBAR.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,t;
  int interval=1;
  int d1=0,d2=0,num;
  double f1,f2;
  double *pmf,pmf_min;

  double din;

  int n_sim;
  int *n;

  double *fene;
  double ***enek;

  double criteria_BAR=1.0e-7;
  int MAXITE=1000;

  int n_total=0;
  double *W,*WT;

  int *frame;
  double *max,*min;
  double *width;
  double ***twod_data,*hist;

  double *c1,*c2;

  char *progname;
  char *enekfilename,*eneklistfilename,*fenefilename;
  char *datalistfilename,*datafilename,*pmffilename;

  FILE *enekfile,*eneklistfile,*fenefile;
  FILE *datalistfile,*datafile,*pmffile;

  int ncid,fene_varid;
  size_t start[1],count[1];

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"hp:i:t:x:"))!=-1) {
    switch(c) {
    case 'i':
      interval=atoi(optarg);
      break;
    case 't':
      criteria_BAR=atof(optarg);
      break;
    case 'x':
      MAXITE=atoi(optarg);
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

  frame=(int *)gcemalloc(sizeof(int)*2);
  max=(double *)gcemalloc(sizeof(double)*2);
  min=(double *)gcemalloc(sizeof(double)*2);
  width=(double *)gcemalloc(sizeof(double)*2);

  argc-=optind;
  argv+=optind;

  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }
  n_sim = atoi(*argv);
  width[0] = atof(*++argv);
  width[1] = atof(*++argv);
  fenefilename = *++argv;
  eneklistfilename = *++argv;
  datalistfilename = *++argv;
  pmffilename = *++argv;

  fene=(double *)gcemalloc(sizeof(double)*n_sim);
  start[0]=0;
  count[0]=n_sim;
  enc_open(fenefilename,NC_NOWRITE,&ncid);
  nc_inq_varid(ncid,"fene",&fene_varid);
  nc_get_vara_double(ncid,fene_varid,start,count,fene);
  encclose(ncid);

  n=(int *)gcemalloc(sizeof(int)*n_sim);
  enek=(double ***)gcemalloc(sizeof(double **)*n_sim);
  for (i=0;i<n_sim;++i) enek[i]=(double **)gcemalloc(sizeof(double *)*n_sim);

  for (i=0;i<n_sim;++i) for (j=0;j<n_sim;++j) enek[i][j]=(double *)gcemalloc(sizeof(double));

  twod_data=(double ***)gcemalloc(sizeof(double **)*2);
  for (i=0;i<2;++i)
    twod_data[i]=(double **)gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<2;++i) for (j=0;j<n_sim;++j) twod_data[i][j]=(double *)gcemalloc(sizeof(double)*1);

  eneklistfile=efopen(eneklistfilename,"r");
  for (i=0;i<n_sim;++i) {
    for (j=0;j<n_sim;++j) {
      getline(&line,&len,eneklistfile);
      line[strlen(line)-1]='\0';
      enekfile=efopen(line,"r");
      getline(&line,&len,enekfile);
      k=0;
      num=0;
      d1 = 1;
      while ( d1 != -1  )  {
	d1=fscanf(enekfile,"%lf",&f1);
	d1=fscanf(enekfile,"%lf",&f1);
	if (k%interval == 0) {
	  enek[j][i]=(double *)gcerealloc(enek[j][i],sizeof(double)*(num+1));
	  enek[j][i][num]=f1;
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
    d1 = 1;
    d2 = 1;
    while ( d1 != -1 && d2 != -1 )  {
      d1=fscanf(datafile,"%lf",&f1);
      d2=fscanf(datafile,"%lf",&f2);
      if (k%interval == 0) {
	twod_data[0][i]=(double *)gcerealloc(twod_data[0][i],sizeof(double)*(num+1));
	twod_data[1][i]=(double *)gcerealloc(twod_data[1][i],sizeof(double)*(num+1));
	twod_data[0][i][num]=f1;
	twod_data[1][i][num]=f2;
	++num;
      }
      ++k;
    }
    fclose(datafile);
  }
  fclose(datalistfile);

  n_total=0;for (i=0;i<n_sim;++i) n_total+=n[i];

  hist=MBAR_AVE_twod(fene,enek,n_sim,n,twod_data,width,max,min,frame);

  pmf=(double *)gcemalloc(sizeof(double)*frame[0]*frame[1]);

  for (i=0;i</*=*/frame[0];++i) {
    for (j=0;j</*=*/frame[1];++j) { 
      if (hist[i*frame[1]+j]!=0.0) {
	pmf[i*frame[1]+j]=/*-*/log(hist[i*frame[1]+j]);
      }
    }
  }

  pmf_min=pmf[0];
  for (i=0;i</*=*/frame[0];++i) {
    for (j=0;j</*=*/frame[1];++j) {
      if ( pmf_min < pmf[i*frame[1]+j] && pmf[i*frame[1]+j]!=0 )
	pmf_min=pmf[i*frame[1]+j];
    }
  }

  for (i=0;i</*=*/frame[0];++i) {
    for (j=0;j</*=*/frame[1];++j) {
      if (hist[i*frame[1]+j]!=0.0) {
  	pmf[i*frame[1]+j]-=pmf_min;
      }
    }
  }

  pmffile=efopen(pmffilename,"w");
  for (i=0;i<frame[0];++i) {
    for (j=0;j<frame[1];++j) {
      if (pmf[i*frame[1]+j]!=0.0)
	fprintf(pmffile,"%10.4lf %10.4lf %10.4lf\n",min[0]+width[0]*i,min[1]+width[1]*j,pmf[i*frame[1]+j]);
      else {
	fprintf(pmffile,"%10.4lf %10.4lf 0.0\n",min[0]+width[0]*i,min[1]+width[1]*j);
      }
    }
  }
  fclose(pmffile);





  return 0;
}

void USAGE(char *progname) {
  printf("[-m fenefilename] -- w mbar ite\n");
  printf("[-t criteria_BAR ] -- criteria of MBAR iteration (default 10-7)\n");
  printf("[-x max_ite_step_of_BAR ] -- max num of MBAR iteration\n");
  printf("[-p] -- interval\n");
  printf("[-h] -- help\n");
  printf("USAGE: %s  n_sim widthx widthy eneklistfilename datalistfilename pmffilename \n", progname);
}

