
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <math.h>
#include <ctype.h>
#include <netcdf.h>
#include <getopt.h>

#include "Optimiza_BFGS_SBRW.h"
#include "Simpson_integ_SBRW.h"
#include "EF.h"

#define ON  0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,ii,jj,num,numnum,dummy;
  double f,dx;

  double k_B=1.98723e-3;             
  double UNITT=418.4070;
  double T=300;

  int interval=1;

  struct data dat;

  int n_sum;   
  double minx, maxx;

  double *w_k,*g_k;
  
  int N_bin=400;
  double x,prob,p;

  double S,sum;

  int periodicflag=ON;
  double delta;
  double periodicity,minv,maxv;

  double *bk_i;

  int num_Simpson=100000;

  double pi;

  char *inputfilelistname,*metadatafilename;
  char *outputfilename,*pmffilename;

  FILE *inputfile,*inputfilelist,*metadatafile;
  FILE *outputfile,*pmffile;
  
  char *line;
  size_t len=0;
  
  int c,d;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;
  
  pi=acos(-1.0);
  
  struct option long_opt[] = {
    {"interval",1,NULL,'m'},
    {"bin",1,NULL,'b'},
    {"num_Simpson",1,NULL,'n'},
    {"T",1,NULL,'t'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"hn:m:t:b:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'b':
      N_bin=atoi(optarg); break;
    case 'n':
      num_Simpson=atoi(optarg);  break;
    case 't':
      T=atof(optarg);     break;
    case 'm':
      interval=atoi(optarg);     break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);  exit(1);
    }
  }

  progname=*argv;  argc-=optind;  argv+=optind;

  if (argc < 7) {
    USAGE(progname);
    exit(1);
  }
  dat.K             = atoi(*argv);
  dat.h             = atof(*++argv);
  dat.n_sim         = atoi(*++argv);
  inputfilelistname = *++argv;
  metadatafilename  = *++argv;
  outputfilename    = *++argv;
  pmffilename       = *++argv;

  periodicflag=ON; periodicity=2.0*pi;
  minv=-1.0*pi; maxv=1.0*pi;

  dat.beta=1.0/(k_B*T/**UNITT*/);

  dat.x_ij=(double **)gcemalloc(sizeof(double *)*(dat.n_sim));
  for (i=0;i<dat.n_sim;++i) dat.x_ij[i]=(double *)gcemalloc(sizeof(double));
  dat.n=(int *)gcemalloc(sizeof(int)*(dat.n_sim));

  n_sum=0;
  inputfilelist=efopen(inputfilelistname,"r");
  for (i=0;i<dat.n_sim;++i) {
    getline(&line,&len,inputfilelist);
    l=strlen(line);
    line[l-1]='\0';

    inputfile=efopen(line,"r");
    num=0; numnum=0;
    d = 1;
    while ( d != -1  )  {
      fscanf(inputfile,"%d",&dummy);
      d=fscanf(inputfile,"%lf",&f);
      if (numnum%interval==0) {
	dat.x_ij[i]=(double *)gcerealloc(dat.x_ij[i],sizeof(double)*(num+1));
	if (periodicflag==ON) {
	  while (f<minv) {
	    f+=periodicity;
	  }
	  while (f>maxv) {
	    f-=periodicity;
	  }
	}
	dat.x_ij[i][num]=f;
	++num;
      }
      ++numnum;
    }
    fclose(inputfile);
    dat.n[i]=num-1;
    n_sum+=dat.n[i];
  }
  fclose(inputfilelist);

  dat.k=(double *)gcemalloc(sizeof(double)*dat.n_sim);
  dat.x_i=(double *)gcemalloc(sizeof(double)*dat.n_sim);
  metadatafile=efopen(metadatafilename,"r");
  for (i=0;i<dat.n_sim;++i) {
    fscanf(metadatafile,"%lf",&(dat.x_i[i])); fscanf(metadatafile,"%lf",&(dat.k[i])); 
    if (periodicflag==ON) {
      while (dat.x_i[i]<minv) {
	dat.x_i[i]+=periodicity;
      }
      while (dat.x_i[i]>maxv) {
	dat.x_i[i]-=periodicity;
      }
    }
  }
  fclose(metadatafile);
  
  minx = dat.x_ij[0][0];
  maxx = dat.x_ij[0][0];
  for (i=0;i<dat.n_sim;++i) {
    for (j=0;j<dat.n[i];++j) {
      if (minx > dat.x_ij[i][j]) minx=dat.x_ij[i][j];
      if (maxx < dat.x_ij[i][j]) maxx=dat.x_ij[i][j];
    }
  }

  dat.x_k=(double *)gcemalloc(sizeof(double)*dat.K);
  dx=(maxx-minx)/(dat.K);
  for (i=0;i<(dat.K);++i) {
    dat.x_k[i]=minx+dx*i;
  }
  
  w_k=(double *)gcemalloc(sizeof(double)*dat.K);
  g_k=(double *)gcemalloc(sizeof(double)*dat.K);

  dat.C_ik=(double **)gcemalloc(sizeof(double *)*dat.n_sim);
  for (i=0;i<dat.n_sim;++i) dat.C_ik[i]=(double *)gcemalloc(sizeof(double)*dat.K);

  dat.S_k_ij=(double ***)gcemalloc(sizeof(double **)*dat.K);
  for (i=0;i<dat.K;++i) dat.S_k_ij[i]=(double **)gcemalloc(sizeof(double *)*dat.n_sim);
  for (i=0;i<dat.K;++i) for (j=0;j<dat.n_sim;++j) dat.S_k_ij[i][j]=(double *)gcemalloc(sizeof(double)*dat.n[j]);

  bk_i=(double *)gcemalloc(sizeof(double)*dat.n_sim);
  for (i=0;i<dat.n_sim;++i) bk_i[i]=dat.beta*dat.k[i];

  for (i=0;i<dat.n_sim;++i) {
    for (j=0;j<dat.K;++j) {
      dat.C_ik[i][j]=Simpson_integ_oneD_SBRW_C_SB(num_Simpson,minv,maxv,dat.x_i[i],dat.x_k[j],bk_i[i],dat.h);
    }
  }

  dat.A=(double *)gcemalloc(sizeof(double)*dat.K);

  for (i=0;i<dat.K;++i) {
    dat.A[i]=Simpson_integ_oneD_SBRW_SB(num_Simpson,minv,maxv,dat.x_k[i],dat.h);
    dat.A[i]=1.0/dat.A[i];
  }

  for (k=0;k<dat.K;++k) {
    for (i=0;i<dat.n_sim;++i) {
      for (j=0;j<dat.n[i];++j) {
	delta=fabs(dat.x_ij[i][j]-dat.x_k[k]);
	if (delta>pi) delta=2.0*pi-delta;

	dat.S_k_ij[k][i][j]=exp(-0.5*delta*delta/dat.h);
      }
    }
  }

  optimize_lnL_BFGS_2(g_k,dat);

  for (i=0;i<dat.K;++i) w_k[i]=exp(g_k[i]);
  sum=0.0; for (i=0;i<dat.K;++i) sum+=w_k[i];
  for (i=0;i<dat.K;++i) w_k[i]=w_k[i]/sum;

  sum=0.0;
  for (i=0;i<dat.K;++i) {
    sum+=w_k[i]*dat.A[i]*Simpson_integ_oneD_SBRW_GB(num_Simpson,minv,maxv,dat.x_k[i],dat.h);
  }

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<dat.K;++i) fprintf(outputfile,"  w[%2d] = %e, \n",i+1, w_k[i]);
  fclose(outputfile);

  pmffile=efopen(pmffilename,"w");
  for (i=0;i<N_bin;++i) {
    x=(maxx-minx)/(double)N_bin*i+minx;
    prob=0.0;
    for (j=0;j<dat.K;++j) {
      delta=fabs(x-dat.x_k[j]);
      if (fabs(delta)>pi) delta=2.0*pi-delta;
      //      p=w_k[j]*dat.A[j]*exp(-0.5*delta*delta/dat.h);
      //      prob+=p;
      prob+=w_k[j]*dat.A[j]*exp(-0.5*delta*delta/dat.h);
    }
    fprintf(pmffile,"%10.8lf %10.8lf %10.8lf\n",x,prob,-1.0*dat.beta*log(prob));
  }
  fclose(pmffile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[--iinterval] interval of data input\n");
  printf("[-h] help \n");
  printf("%s [-h] K h n_sim inputfilelistname metadatafilename pmffilename \n",progname);
}
