
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <math.h>
#include <ctype.h>
#include <netcdf.h>
#include <getopt.h>

#include "Optimiza_BFGS_2.h"
#include "EF.h"

#define ON  0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,ii,jj,k,l,num,numnum,dummy;
  double f,dx;

  double k_B=1.98723e-3;             
  double UNITT=418.4070;
  double T=300;

  int interval=1;

  struct data dat;

  int n_sum;   
  double minx, maxx;

  double *f_i,*g_i;
  
  int N_bin=400;
  double x,prob;

  double S,sum;

  double *bk_i;

  int periodicflag=OFF;
  double delta;
  double periodicity,minv,maxv;

  double pi;

  char *inputfilelistname,*metadatafilename;
  char *pmffilename;

  FILE *inputfile,*inputfilelist,*metadatafile;
  FILE *pmffile;
  
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
    {"pi",0,NULL,'i'},
    {"T",1,NULL,'t'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"him:t:b:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'b':
      N_bin=atoi(optarg); break;
    case 'i':
      periodicflag=ON; periodicity=2.0*pi;
      minv=-1.0*pi; maxv=1.0*pi;     break;
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  dat.n_sim         = atoi(*argv);
  inputfilelistname = *++argv;
  metadatafilename  = *++argv;
  pmffilename       = *++argv;

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
  dat.K=dat.n_sim;

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
  
  f_i=(double *)gcemalloc(sizeof(double)*dat.K);
  g_i=(double *)gcemalloc(sizeof(double)*dat.K);

  dat.exp_kii_xij=(double ***)gcemalloc(sizeof(double **)*dat.n_sim);
  for (i=0;i<dat.n_sim;++i) 
    dat.exp_kii_xij[i]=(double **)gcemalloc(sizeof(double *)*dat.n_sim);
  for (i=0;i<dat.n_sim;++i) 
    for (j=0;j<dat.n_sim;++j) 
      dat.exp_kii_xij[i][j]=(double *)gcemalloc(sizeof(double)*dat.n[j]);

  bk_i=(double *)gcemalloc(sizeof(double)*dat.n_sim);
  for (i=0;i<dat.n_sim;++i) bk_i[i]=dat.beta*dat.k[i];

  for (i=0;i<dat.n_sim;++i){
    for (ii=0;ii<dat.n_sim;++ii){
      for (jj=0;jj<dat.n[ii];++jj){
	delta=fabs(dat.x_ij[ii][jj]-dat.x_i[i]);
	/*********************************************/
        /* if (i==17 && ii==17 && jj==0) {	     */
	/*   printf("%2d-%2d-%2d : d=%3lf\n",	     */
	/* 	 i,ii,jj,delta);		     */
	/* }					     */
        /*********************************************/
	if (periodicflag==ON) {
	  if (fabs(delta)>0.5*periodicity) {
	    delta=periodicity-delta;
	  }
	}
	dat.exp_kii_xij[i][ii][jj]=exp(-1.0*(0.5*bk_i[i]*delta*delta));
	/**********************************************************************************/
        /* if (i==17 && ii==17 && jj==0) {						  */
	/*   printf("%2d-%2d-%2d : c=%3lf d=%3lf d*d=%3lf b=%3lf\n",			  */
	/* 	 i,ii,jj,dat.exp_kii_xij[i][ii][jj],delta,delta*delta,bk_i[i]);		  */
	/* }										  */
        /**********************************************************************************/
      }
    }
  }

  optimize_lnL_BFGS(g_i,dat);

  for (i=0;i<dat.n_sim;++i) f_i[i]=exp(g_i[i]);

  pmffile=efopen(pmffilename,"w");
  l=0;
  for (i=0;i<dat.n_sim;++i) {
    for (j=0;j<dat.n[i];++j) {      
      prob=0.0;
      for (ii=0;ii<dat.n_sim;++ii) prob+=dat.n[ii]*f_i[ii]*dat.exp_kii_xij[ii][i][j];
      prob=1.0/prob;

      fprintf(pmffile,"%10.8lf %10.8lf %10.8lf\n",dat.x_ij[i][j],prob,-1.0*log(prob));
      ++l;
    }
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
