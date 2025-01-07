
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "IO.h"

#include "MBAR.h"

#define ON 1
#define OFF 0

#define UP 1
#define DOWN 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,t;
  int flag=OFF,iflag=OFF,jflag=OFF;
  int interval=1,ueneinterval=1;
  int d=0,num;
  double sum=0.0;
  double f,f1,f2;

  int MODE=UP;

  double din;

  int n_sim,numk=0;
  int *n;

  double T=300,*T_sim;
  double k_B=1.98723e-3;
  double beta;

  double *fene;
  double ***enek,**expeneraw;

  int framex,framey;
  double maxx,minx,maxy,miny;
  double widthx,widthy;
  double ***td_data,**hist,**pmf;
  double pmf_min;

  double criteria_BAR=1.0e-7;
  int MAXITE=1000;

  char *progname;
  char *enekfilename,*eneklistfilename,*fenefilename;
  char *datalistfilename,*datafilename;
  char *histfilename,*pmffilename;
  char *Tsimlistfilename;

  FILE *enekfile,*eneklistfile,*fenefile;
  FILE *datalistfile,*datafile;
  FILE *histfile,*pmffile;
  FILE *Tsimlistfile;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"hdijmt:p:u:"))!=-1) {
    switch(c) {
    case 'm':
      flag=ON;
      break;
    case 'd':
      MODE=DOWN;
      break;
    case 'i':
      iflag=ON;
      break;
    case 'j':
      jflag=ON;
      break;
    case 'p':
      interval=atoi(optarg);
      break;
    case 'u':
      ueneinterval=atoi(optarg);
      break;
    case 't':
      T=atof(optarg);
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

  if (argc < 10) {
    USAGE(progname);
    exit(1);
  }
  numk = atoi(*argv)-1;
  n_sim = atoi(*++argv);
  widthx = atof(*++argv);
  widthy = atof(*++argv);
  fenefilename = *++argv;
  eneklistfilename = *++argv;
  Tsimlistfilename = *++argv;
  datalistfilename = *++argv;
  histfilename = *++argv;
  pmffilename = *++argv;

  fene=(double *)gcemalloc(sizeof(double)*n_sim);
  if (flag==OFF) {
    fenefile=efopen(fenefilename,"r");
    getline(&line,&len,fenefile);
    for (i=0;i<n_sim;++i) {
      fscanf(fenefile,"%d %lf",&d,&fene[i]);
      fene[i]=exp(fene[i]);
    }
    fclose(fenefile);
  }

  n=(int *)gcemalloc(sizeof(int)*n_sim);
  enek=(double ***)gcemalloc(sizeof(double **)*n_sim);
  for (i=0;i<n_sim;++i) {
    enek[i]=(double **)gcemalloc(sizeof(double *)*n_sim);
  }
  for (i=0;i<n_sim;++i) {
    for (j=0;j<n_sim;++j) {
      enek[i][j]=(double *)gcemalloc(sizeof(double ));
    }
  }
  td_data=(double ***)gcemalloc(sizeof(double **)*n_sim);
  for (i=0;i<n_sim;++i) td_data[i]=(double **)gcemalloc(sizeof(double *));
  expeneraw=(double **)gcemalloc(sizeof(double *)*n_sim);
  T_sim=(double *)gcemalloc(sizeof(double)*n_sim);

  Tsimlistfile=efopen(Tsimlistfilename,"r");
  for (i=0;i<n_sim;++i) {
    fscanf(Tsimlistfile,"%lf",&T_sim[i]);
  }
  fclose(Tsimlistfile);

  eneklistfile=efopen(eneklistfilename,"r");
  for (i=0;i<n_sim;++i) {
    getline(&line,&len,eneklistfile);
    l=strlen(line);
    line[l-1]='\0';
    enekfile=efopen(line,"r");
    getline(&line,&len,enekfile);
    k=0;
    num=0;
    d = 1;
    while ( d != -1  )  {
      d=fscanf(enekfile,"%d",&l);
      d=fscanf(enekfile,"%lf",&f);
      if (k%ueneinterval == 0) {
	expeneraw[i]=(double *)gcerealloc(expeneraw[i],sizeof(double)*(/*k*/num+1));
	expeneraw[i][num]=f;
	++num;
      }
      ++k;
    } 
    fclose(enekfile);
    n[i]=num-1;
  }
  fclose(eneklistfile);

  for (i=0;i<n_sim;++i) {
    beta=1.0/k_B*(1.0/T_sim[i]-1.0/T);
    for (j=0;j<n_sim;++j) {
      num=0;
      for (k=0;k<n[j];++k) {
	f=expeneraw[j][k];
	enek[i][j]=(double *)gcerealloc(enek[i][j],sizeof(double)*(num+1));
	enek[i][j][num]=exp(-1.0*f*beta);
	++num;
      } 
    }
  }

  datalistfile=efopen(datalistfilename,"r");
  for (i=0;i<n_sim;++i) {
    getline(&line,&len,datalistfile);
    line[strlen(line)-1]='\0';
    datafile=efopen(line,"r");
    k=0;
    num=0;
    d = 1;
    while ( d != -1  )  {
      if (jflag==ON)
	d=fscanf(datafile,"%d",&l);
      d=fscanf(datafile,"%lf %lf",&f1,&f2);
      if (k%interval == 0) {
	td_data[i]=(double **)gcerealloc(td_data[i],sizeof(double *)*(num+1));
	td_data[i][num]=(double *)gcemalloc(sizeof(double)*2);
	td_data[i][num][0]=f1;
	td_data[i][num][1]=f2;
	++num;
      }
      ++k;
    } 
    fclose(datafile);
  }
  fclose(datalistfile);

  if (flag==ON) MBAR_ite_high_speed_2(fene,enek,n_sim,n,criteria_BAR,MAXITE);
    
  hist=MBAR_AVE_twod_multi_temp(fene,enek,n_sim,n,td_data,widthx,widthy,&maxx,&minx,&maxy,&miny,&framex,&framey,numk);

  pmf=(double **)gcemalloc(sizeof(double *)*framex);
  for (i=0;i<=framex;++i) pmf[i]=(double *)gcemalloc(sizeof(double)*framey);

  pmf_min=0.0;
  for (i=0;i</*=*/framex;++i) {
    for (j=0;j</*=*/framey;++j) {
      if ( pmf_min == 0.0 && hist[i][j] > 0.0 ) pmf_min=hist[i][j];
      if ( MODE == UP ) {
	if ( pmf_min < hist[i][j] && hist[i][j] > 0.0 ) pmf_min=hist[i][j];
      }
      else {
	if ( pmf_min > hist[i][j] && hist[i][j] > 0.0 ) pmf_min=hist[i][j];
      }
    }
  }

  for (i=0;i</*=*/framex;++i) {
    for (j=0;j</*=*/framey;++j) { 
      if (hist[i][j]!=0.0) {
      if ( MODE == UP ) {
	pmf[i][j]=-log(hist[i][j])+log(pmf_min);
      }
      else {
	pmf[i][j]=-(/*-*/log(hist[i][j])/*+*/-log(pmf_min));
      }  
      }
    }
  }

  if (flag==ON) {
    fenefile=efopen(fenefilename,"w");
    fprintf(fenefile,"# fene\n");
    for (i=0;i<n_sim;++i)
      fprintf(fenefile,"%d %12.10lf\n",i,fene[i]);
    fclose(fenefile);
  }

  histfile=efopen(histfilename,"w");
  fprintf(histfile,"x y hist \n");
  for (i=0;i<framex;++i)
    for (j=0;j<framey;++j)
      fprintf(histfile,"%10.4lf %10.4lf %10.4lf\n",minx+widthx*i,miny+widthy*i,hist[i][j]);
  fclose(histfile);

  pmffile=efopen(pmffilename,"w");
  //  fprintf(pmffile,"x y  pmf  \n");
  for (i=0;i<framex;++i) {
    for (j=0;j<framey;++j) {
      if (pmf[i][j]!=0.0)
	fprintf(histfile,"%10.4lf %10.4lf %10.4lf\n",minx+widthx*i,miny+widthy*j,pmf[i][j]);
      else {
	if ( MODE == UP )
	  fprintf(histfile,"%10.4lf %10.4lf Nan\n",minx+widthx*i,miny+widthy*j);
	else
	  fprintf(histfile,"%10.4lf %10.4lf 0.0\n",minx+widthx*i,miny+widthy*j);
      }
    }
  }
  fclose(pmffile);

  printf("%d %d\n",framex,framey);
  
  return 0;
}

void USAGE(char *progname) {
  printf("-p -- interval\n");
  printf("-h -- help\n");
  printf("USAGE: %s n_sim width fenefilename enelistfilename datalistfilename histfilename histfilename pmffilename errorfile perrorfilename \n", progname);
}

