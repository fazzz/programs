
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "netcdf_mine.h"
#include "EF.h"
#include "IO.h"

#include "MBAR.h"

#define ON 1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int bflag=OFF;
  int interval=1;
  int d=0,num;
  int num_sim;
  int *numstep;
  double f;
  double ***ene;
  double *fene;

  double criteria_BAR=1.0e-7;
  int MAXITE=1000;

  char *progname;
  char *enefilename,*enelistfilename,*fenefilename;
  FILE *enefile,*enelistfile,*fenefile;

  int ncid,indexofsim_dimid,fene_varid;
  size_t start[1],count[1];

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"hnp:c:x:"))!=-1) {
    switch(c) {
    case 'n':
      bflag=ON;
      break;
    case 'c':
      criteria_BAR=atof(optopt);
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

  argc-=optind;
  argv+=optind;

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  num_sim = atoi(*argv);
  enelistfilename  = *++argv;
  fenefilename = *++argv;

  numstep=(int *)gcemalloc(sizeof(int)*num_sim);
  ene=(double ***)gcemalloc(sizeof(double **)*num_sim);
  for (i=0;i<num_sim;++i) {
    ene[i]=(double **)gcemalloc(sizeof(double *)*num_sim);
  }
  for (i=0;i<num_sim;++i) {
    for (j=0;j<num_sim;++j) {
      ene[i][j]=(double *)gcemalloc(sizeof(double));
    }
  }

  enelistfile=efopen(enelistfilename,"r");
  for (i=0;i<num_sim;++i) {
    for (j=0;j<num_sim;++j) {
      getline(&line,&len,enelistfile);
      l=strlen(line);
      line[l-1]='\0';
      enefile=efopen(line,"r");
      getline(&line,&len,enefile);
      k=0;
      num=0;
      d = 1;
      while ( d != -1  )  {
	d=fscanf(enefile,"%d",&l);
	d=fscanf(enefile,"%lf",&f);
	if (k%interval == 0) {
	  ene[j][i]=(double *)gcerealloc(ene[j][i],sizeof(double)*(k+1));
	  ene[j][i][k]=f;
	  ++num;
	}
	++k;
      } 
      fclose(enefile);
      numstep[i]=num-1;
    }
  }
  fclose(enelistfile);

  fene=(double *)gcemalloc(sizeof(double)*num_sim);
  
  MBAR_ite(fene,ene,num_sim,numstep,criteria_BAR,MBAR_ite);

  if (bflag==OFF) {
    fenefile=efopen(fenefilename,"w");
    fprintf(fenefile,"# fene \n");
    for (i=0;i<num_sim;++i)
      fprintf(fenefile,"%3d %10.8lf\n",i,fene[i]);
    fclose(fenefile);
  }
  else {
    enc_create(fenefilename,NC_CLOBBER,&ncid);
    enc_def_dim(ncid,"indexofsim",NC_UNLIMITED,&indexofsim_dimid);
    enc_def_var(ncid,"fene",NC_DOUBLE,1,indexofsim_dimid,&fene_varid);
    nc_put_att_text(ncid,fene_varid,UNITS,strlen(ENE_UNIT),ENE_UNIT);
    nc_enddef(ncid);
    
    start[0]=0;
    count[0]=num_sim;
    nc_put_vara_double(ncid,fene_varid,&start,&count,fene);
    encclose(ncid);
  }

  return 0;
}

void USAGE(char *progname) {
  printf("-n -- netcdf output  \n");
  printf("-c -- criteria of iterarion of BAR (default 10^-07)  \n");
  printf("-x -- max num of iteration (default 1000)\n");
  printf("-h -- help\n");
  printf("USAGE: [-n] [-c] [-x] [-h] %s num_sim enelistfilename fenefilename \n", progname);
}

