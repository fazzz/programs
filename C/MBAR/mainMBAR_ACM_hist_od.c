
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

#define ON 1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,t;
  int aflag=OFF,bflag=OFF,mflag=OFF;
  int interval=1;
  int d=0,num;
  double f;

  double din;

  int n_sim;
  int *n;

  double *fene;
  double ***enek;

  double criteria_BAR=1.0e-7;
  int MAXITE=1000;

  int n_total=0;
  double *W,*WT,*W2,*WT2;

  int frame;
  double max,min;
  double width;
  double **od_data,*hist,*hist_error,*hist_error2,*covm;

  double *c1,*c2;

  char *progname;
  char *enekfilename,*eneklistfilename,*fenefilename;
  char *datalistfilename,*datafilename;
  char *histfilename,*histfilename2,*ehistfilename;

  FILE *enekfile,*eneklistfile,*fenefile;
  FILE *datalistfile,*datafile;
  FILE *histfile,*histfile2,*ehistfile,*logfile;

  int ncid,fene_varid;
  size_t start[1],count[1];

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"hap:n:m:t:x:"))!=-1) {
    switch(c) {
    case 'a':
      aflag=ON;
      break;
    case 'n':
      bflag=ON;
      fenefilename=optarg;
      break;
    case 'm':
      mflag=ON;
      fenefilename=optarg;
      break;
    case 'p':
      interval=atoi(optarg);
      break;
    case 't':
      criteria_BAR=atof(optarg);
      break;
    case 'x':
      MAXITE=atoi(optarg);
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  if ( bflag==ON && mflag==ON ) {
    printf("option error: -b -m");
    exit(1);
  }

  argc-=optind;
  argv+=optind;

  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }
  n_sim = atoi(*argv);
  width = atof(*++argv);
  eneklistfilename = *++argv;
  datalistfilename = *++argv;
  histfilename = *++argv;
  ehistfilename = *++argv;

  fene=(double *)gcemalloc(sizeof(double)*n_sim);
  if (bflag==OFF && mflag==OFF) {
    fenefile=efopen(fenefilename,"r");
    getline(&line,&len,fenefile);
    for (i=0;i<n_sim;++i)
      fscanf(fenefile,"%d %lf",&d,&fene[i]);
    fclose(fenefile);
  }
  if (bflag==ON) {
    start[0]=0;
    count[0]=n_sim;
    enc_open(fenefilename,NC_NOWRITE,&ncid);
    nc_inq_varid(ncid,"fene",&fene_varid);
    nc_get_vara_double(ncid,fene_varid,start,count,fene);
    encclose(ncid);
   }

  n=(int *)gcemalloc(sizeof(int)*n_sim);
  enek=(double ***)gcemalloc(sizeof(double **)*n_sim);
  for (i=0;i<n_sim;++i) enek[i]=(double **)gcemalloc(sizeof(double *)*n_sim);

  for (i=0;i<n_sim;++i) for (j=0;j<n_sim;++j) enek[i][j]=(double *)gcemalloc(sizeof(double));

  od_data=(double **)gcemalloc(sizeof(double *)*n_sim);
  for (i=0;i<n_sim;++i) od_data[i]=(double **)gcemalloc(sizeof(double));

  eneklistfile=efopen(eneklistfilename,"r");
  for (i=0;i<n_sim;++i) {
    for (j=0;j<n_sim;++j) {
      getline(&line,&len,eneklistfile);
      line[strlen(line)-1]='\0';
      enekfile=efopen(line,"r");
      getline(&line,&len,enekfile);
      k=0;
      num=0;
      d = 1;
      while ( d != -1  )  {
	d=fscanf(enekfile,"%lf",&f);
	d=fscanf(enekfile,"%lf",&f);
	if (k%interval == 0) {
	  enek[j][i]=(double *)gcerealloc(enek[j][i],sizeof(double)*(num+1));
	  enek[j][i][num]=f;
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
    d = 1;
    while ( d != -1  )  {
      d=fscanf(datafile,"%lf",&f);
      if (k%interval == 0) {
	od_data[i]=(double *)gcerealloc(od_data[i],sizeof(double)*(num+1));
	od_data[i][num]=f;
	++num;
      }
      ++k;
    } 
    fclose(datafile);
  }
  fclose(datalistfile);

  if (mflag==ON) MBAR_ite(fene,enek,n_sim,n,criteria_BAR,MAXITE);

  n_total=0;for (i=0;i<n_sim;++i) n_total+=n[i];

  hist=MBAR_AVE_oned(fene,enek,n_sim,n,od_data,width,&max,&min,&frame);
  histfile=efopen(histfilename,"w");
  fprintf(histfile,"# hist \n");
  for (i=0;i<frame;++i)
    fprintf(histfile,"%10.4e %10.4e\n",min+width*i,hist[i]);
  fclose(histfile);

  if (aflag==OFF) {
    W=gcemalloc(sizeof(double)*n_total*(n_sim+2)); 
    WT=(double *)gcemalloc(sizeof(double)*n_total*(n_sim+2));

    l=-1;
    for (i=0;i<n_sim;++i) {
      for (t=0;t<n[i];++t) {
	++l;
	for (j=0;j<n_sim;++j) {
	  din=0.0;
	  for (k=0;k<n_sim;++k) {
	    din+=n[k]*exp(fene[k]-enek[k][i][t]);
	  }
	  W[l*(n_sim+2)+j]=exp(fene[j]-enek[j][i][t])/din;
	  WT[j*n_total+l] =W[l*(n_sim+2)+j];
	}
      }
    }

    hist_error=MBAR_ACM2_histod(fene,enek,n_sim,n,od_data,width,&max,&min,&frame,W,WT);

    ehistfile=efopen(ehistfilename,"w");
    fprintf(ehistfile,"# ehist \n");
    for (i=0;i<frame;++i)
      fprintf(ehistfile,"%10.4e %10.4e\n",width*i,hist_error[i]);
    fclose(ehistfile);

    W2=gcemalloc(sizeof(double)*n_total*(n_sim)); 
    WT2=(double *)gcemalloc(sizeof(double)*n_total*(n_sim));

    l=-1;
    for (i=0;i<n_sim;++i) {
      for (t=0;t<n[i];++t) {
	++l;
	for (j=0;j<n_sim;++j) {
	  din=0.0;
	  for (k=0;k<n_sim;++k) {
	    din+=n[k]*exp(fene[k]-enek[k][i][t]);
	  }
	  W2[l*n_sim+j]=exp(fene[j]-enek[j][i][t])/din;
	  WT2[j*n_total+l] =W2[l*n_sim+j];
	}
      }
    }

    covm=(double *)gcemalloc(sizeof(double)*n_sim*n_sim);
    covm=MBAR_ACM2(n_sim,n,W2,WT2);

    logfile=efopen("logMBAR_error_est.txt","w");
    fprintf(logfile,"# estimated error \n");
    for (i=0;i<n_sim;++i)
      fprintf(logfile,"%d %10.4lf\n",i,covm[i*n_sim+i]-2.0*covm[i*n_sim]+covm[0]);
    fclose(logfile);
  }  

  return 0;
}

void USAGE(char *progname) {
  printf("[-a ] -- wo error calculation\n");
  printf("[-n fenefilename] -- nc input\n");
  printf("[-m fenefilename] -- w mbar ite\n");
  printf("[-t criteria_BAR ] -- criteria of MBAR iteration (default 10-7)\n");
  printf("[-x max_ite_step_of_BAR ] -- max num of MBAR iteration\n");
  printf("[-p] -- interval\n");
  printf("[-h] -- help\n");
  printf("USAGE: %s [-a] [-n fenefilename] [-m fenefilename] [-p] [-h] n_sim width eneklistfilename datalistfilename histfilename ehistfilename\n", progname);
}
