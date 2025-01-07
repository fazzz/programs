
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <math.h>
#include <ctype.h>
#include <netcdf.h>
#include <getopt.h>

#include "Optimize_BFGS_TPBRW.h"
#include "Simpson_integ_TPBRW.h"
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

  double *g_k;
  double a_0,*a_k,*b_k;
  
  int N_bin=400;
  double x,free_ene,p;

  double S,sum;

  double delta;
  double periodicity,/*minv,maxv,*/periodicflag=ON;

  //  double *bk_i;

  double k_x;

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

  (dat).num_Simpson=100000;
  
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
      dat.num_Simpson=atoi(optarg);  break;
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

  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }
  dat.K             = atoi(*argv);
  dat.n_sim         = atoi(*++argv);
  inputfilelistname = *++argv;
  metadatafilename  = *++argv;
  outputfilename    = *++argv;
  pmffilename       = *++argv;

  periodicity=2.0*pi;
  dat.minv=-1.0*pi; 
  dat.maxv=1.0*pi;

  dat.beta=1.0/(k_B*T/*UNIT*/);

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
	  while (f<dat.minv) {
	    f+=periodicity;
	  }
	  while (f>dat.maxv) {
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
      while (dat.x_i[i]<dat.minv) {
	dat.x_i[i]+=periodicity;
      }
      while (dat.x_i[i]>dat.maxv) {
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

  /******************************************************/
  /* dat.x_k=(double *)gcemalloc(sizeof(double)*dat.K); */
  /* dx=(maxx-minx)/(dat.K);			        */
  /* for (i=0;i<(dat.K);++i) {			        */
  /*   dat.x_k[i]=minx+dx*i;			        */
  /* }						        */
  /******************************************************/
  
  /**************************************************/
  /* a_k=(double *)gcemalloc(sizeof(double)*dat.K); */
  /* b_k=(double *)gcemalloc(sizeof(double)*dat.K); */
  /**************************************************/

  /**********************************************************/
  /* dat.C_i=(double *)gcemalloc(sizeof(double)*dat.n_sim); */
  /**********************************************************/

  /****************************************************************************************/
  /* dat.C_sin_ik=(double **)gcemalloc(sizeof(double *)*dat.n_sim);			  */
  /* for (i=0;i<dat.n_sim;++i) dat.C_sin_ik[i]=(double *)gcemalloc(sizeof(double)*dat.K); */
  /* 											  */
  /* dat.C_cos_ik=(double **)gcemalloc(sizeof(double *)*dat.n_sim);			  */
  /* for (i=0;i<dat.n_sim;++i) dat.C_cos_ik[i]=(double *)gcemalloc(sizeof(double)*dat.K); */
  /****************************************************************************************/

  dat.sin_k_ij=(double ***)gcemalloc(sizeof(double **)*dat.K);
  for (i=0;i<dat.K;++i) dat.sin_k_ij[i]=(double **)gcemalloc(sizeof(double *)*dat.n_sim);
  for (i=0;i<dat.K;++i) for (j=0;j<dat.n_sim;++j) dat.sin_k_ij[i][j]=(double *)gcemalloc(sizeof(double)*dat.n[j]);

  dat.cos_k_ij=(double ***)gcemalloc(sizeof(double **)*dat.K);
  for (i=0;i<dat.K;++i) dat.cos_k_ij[i]=(double **)gcemalloc(sizeof(double *)*dat.n_sim);
  for (i=0;i<dat.K;++i) for (j=0;j<dat.n_sim;++j) dat.cos_k_ij[i][j]=(double *)gcemalloc(sizeof(double)*dat.n[j]);

  dat.bk_i=(double *)gcemalloc(sizeof(double)*dat.n_sim);
  for (i=0;i<dat.n_sim;++i) dat.bk_i[i]=dat.beta*dat.k[i];

  /**********************************************************************************************************/
  /* for (i=0;i<dat.n_sim;++i) {									     */
  /*     dat.C_i[i]=Simpson_integ_oneD_TPBRW_C(num_Simpson,minv,maxv,dat.x_i[i],bk_i[i]);		     */
  /*   for (j=0;j<dat.K;++j) {										     */
  /*         dat.C_sin_ik[i][j]=Simpson_integ_oneD_TPBRW_C_sin(num_Simpson,minv,maxv,dat.x_i[i],j,bk_i[i]); */
  /*         dat.C_cos_ik[i][j]=Simpson_integ_oneD_TPBRW_C_cos(num_Simpson,minv,maxv,dat.x_i[i],j,bk_i[i]); */
  /*   }												     */
  /* }													     */
  /**********************************************************************************************************/

  /************************************************************************************/
  /* dat.A=(double *)gcemalloc(sizeof(double)*dat.K);				      */
  /* 										      */
  /* for (i=0;i<dat.K;++i) {							      */
  /*   dat.A[i]=Simpson_integ_oneD_TPBRW_TPB(num_Simpson,minv,maxv,dat.x_k[i],dat.h); */
  /*   dat.A[i]=1.0/dat.A[i];							      */
  /* }										      */
  /************************************************************************************/

  for (k=0;k<dat.K;++k) {
    for (i=0;i<dat.n_sim;++i) {
      for (j=0;j<dat.n[i];++j) {
	k_x=(double)(k+1)*dat.x_ij[i][j];
	while (k_x > 2.0*pi) k_x-=2.0*pi;
	while (k_x <= 0.0)   k_x+=2.0*pi;

	dat.sin_k_ij[k][i][j]=sin(k_x);
	dat.cos_k_ij[k][i][j]=cos(k_x);
      }
    }
  }

  g_k=(double *)gcemalloc(sizeof(double)*(dat.K*2+1));
  optimize_lnL_BFGS_TPBRW(g_k,dat);

  /***********************************************/
  /* for (i=0;i<dat.K;++i) w_k[i]=exp(g_k[i]);	 */
  /* sum=0.0; for (i=0;i<dat.K;++i) sum+=w_k[i]; */
  /* for (i=0;i<dat.K;++i) w_k[i]=w_k[i]/sum;	 */
  /***********************************************/

  /************************************************************************************************/
  /* sum=0.0;											  */
  /* for (i=0;i<dat.K;++i) {									  */
  /*   sum+=w_k[i]*dat.A[i]*Simpson_integ_oneD_TPBRW_TPB(num_Simpson,minv,maxv,dat.x_k[i],dat.h); */
  /* }												  */
  /************************************************************************************************/

  a_k=(double *)gcemalloc(sizeof(double)*dat.K);
  b_k=(double *)gcemalloc(sizeof(double)*dat.K);

  a_0=g_k[0];
  for (i=0;i<dat.K;++i) {
    a_k[i]=g_k[i+1];
    b_k[i]=g_k[i+dat.K+1];
  }

  outputfile=efopen(outputfilename,"w");
  fprintf(outputfile,"  a0 = %e, \n", a_0);
  for (i=0;i<dat.K;++i) fprintf(outputfile,"  a%2d = %e, \n",i+1, a_k[i]);
  for (i=0;i<dat.K;++i) fprintf(outputfile,"  b%2d = %e, \n",i+1, b_k[i]);
  fclose(outputfile);

  pmffile=efopen(pmffilename,"w");
  for (i=0;i<N_bin;++i) {
    x=(maxx-minx)/(double)N_bin*i+minx;
    free_ene=0.5*a_0;
    for (j=0;j<dat.K;++j) {
      free_ene+=a_k[j]*sin((j+1)*x)+b_k[j]*cos((j+1)*x);
    }
    fprintf(pmffile,"%10.8lf %10.8lf \n",x,free_ene);
  }
  fclose(pmffile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[--iinterval] interval of data input\n");
  printf("[-h] help \n");
  printf("%s [-h] K n_sim inputfilelistname metadatafilename pmffilename \n",progname);
}
