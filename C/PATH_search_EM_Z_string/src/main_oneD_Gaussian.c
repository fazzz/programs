
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <netcdf.h>
#include <getopt.h>

#include "EM_reweight.h"

#include "Simpson_integ.h"

#include "Gaussian.h"

#include "EF.h"

#define ON 0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,num,dummy;
  double f;

  int num_sim,*n,n_sum;  // # of simulation, # of snapshots
  double **x_ij;         // raw data

  // parameters for umbrella sampling
  double *k_umbrella,*x0;
  
  // parameters of Mixed-Gaussian
  int K; // # of Gaussian 
  double *nyu_k,*Sigma_k, *pi_k;  // mean,variance,weight, # of snapshots for k-th prob.
  double ***gammak_ij, **inv_f_ik, *f_i; // responsibilities
  double L,Lnew,dL=0.0,ddL=1.0,threhold=1.0e-4; // threhold for conversion etc.

  // parameters for Simpson integration
  int num_Simpson=10000; // # of bins (define width of integration)
  double minx, maxx; // min and max values

  int N_bin=50;
  double x,prob;

  double pi;

  char *inputfilelistname,*metadatafilename,*outputfilename,*pmffilename,*pmftempfilename="PMFOUT";
  FILE *inputfile,*inputfilelist,*metadatafile,*outputfile,*pmffile,*pmftempfile;
  
  char *line;
  size_t len=0;
  
  int c,d;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;
  
  pi=acos(-1.0);
  
  struct option long_opt[] = {
    {"minx",1,NULL,'x'},
    {"maxx",1,NULL,'a'},
    {"num_Simpson",1,NULL,'n'},
    {"num_bin",1,NULL,'b'},
    {"K",1,NULL,'K'},
    {"pmftemp",1,NULL,'p'},
    {"tol",1,NULL,'t'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"hb:K:x:a:n:p:t:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'K':
      K = atoi(optarg);  break;
    case 'x':
      minx=atof(optarg);  break;
    case 'a':
      maxx=atof(optarg);  break;
    case 'n':
      num_Simpson=atoi(optarg);  break;
    case 'b':
      N_bin=atoi(optarg);  break;
    case 'p':
      pmftempfilename=optarg;  break;
    case 't':
      threhold=atof(optarg);  break;
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
  num_sim           = atoi(*argv);
  K                 = atoi(*++argv);
  inputfilelistname = *++argv;
  metadatafilename  = *++argv;
  outputfilename    = *++argv;
  pmffilename       = *++argv;

  x_ij=(double **)gcemalloc(sizeof(double *)*(num_sim));
  for (i=0;i<num_sim;++i) x_ij[i]=(double *)gcemalloc(sizeof(double));
  n=(int *)gcemalloc(sizeof(int)*(num_sim));

  n_sum=0;
  inputfilelist=efopen(inputfilelistname,"r");
  for (i=0;i<num_sim;++i) {
    getline(&line,&len,inputfilelist);
    l=strlen(line);
    line[l-1]='\0';

    inputfile=efopen(line,"r");
    num=0;
    d = 1;
    while ( d != -1  )  {
      fscanf(inputfile,"%d",&dummy);
      d=fscanf(inputfile,"%lf",&f);
      x_ij[i]=(double *)gcerealloc(x_ij[i],sizeof(double)*(num+1));
      x_ij[i][num]=f;
      ++num;
    }
    fclose(inputfile);
    n[i]=num-1;
    n_sum+=n[i];
  }
  fclose(inputfilelist);

  k_umbrella=(double *)gcemalloc(sizeof(double)*num_sim);
  x0=(double *)gcemalloc(sizeof(double)*num_sim);
  metadatafile=efopen(metadatafilename,"r");
  for (i=0;i<K;++i) {
    fscanf(metadatafile,"%lf",&x0[i]); fscanf(metadatafile,"%lf",&k_umbrella[i]); 
  }
  fclose(metadatafile);
  
  minx = x_ij[0][0];
  maxx = x_ij[0][0];
  for (i=0;i<num_sim;++i) {
    for (j=0;j<n[i];++j) {
      if (minx > x_ij[i][j]) minx=x_ij[i][j];
      if (maxx < x_ij[i][j]) maxx=x_ij[i][j];
    }
  }

  nyu_k=(double *)gcemalloc(sizeof(double)*K);
  Sigma_k=(double *)gcemalloc(sizeof(double)*K);
  pi_k=(double *)gcemalloc(sizeof(double)*K);

  //Initialization
  for (i=0;i<K;++i) {
    nyu_k[i]=(maxx-minx)/(double)K*i+minx;
    Sigma_k[i]=1.0;
    pi_k[i]=1.0/(double)K;
  }
  gammak_ij=(double ***)gcemalloc(sizeof(double **)*K);
  for (i=0;i<K;++i) {
    gammak_ij[i]=(double **)gcemalloc(sizeof(double *)*num_sim);
    for (j=0;j<num_sim;++j) {
      gammak_ij[i][j]=(double *)gcemalloc(sizeof(double)*n[j]);
    }
  }
  inv_f_ik=(double **)gcemalloc(sizeof(double *)*num_sim);
  for (i=0;i<num_sim;++i) {
    inv_f_ik[i]=(double *)gcemalloc(sizeof(double)*K);
  }
  f_i=(double *)gcemalloc(sizeof(double)*num_sim);

  // main loop. EM algorithum.
  i=0;
  dL=threhold;
  while (ddL >= threhold || i==0) {

    E_step_oneD(num_sim,n,x_ij,              // # of simulation, # of snapshots, data,
		minx,maxx,num_Simpson,       // parameters for Simpson integration 1
		k_umbrella,x0,               // parameters for Simpson integration 2
		K,nyu_k,Sigma_k,pi_k,        // parameters of MG
		gammak_ij,inv_f_ik,f_i,      // responsibilities
		pi);

    M_step_oneD(num_sim, n, n_sum, x_ij,     // # of simulation, # of snapshots, data,
		minx, maxx, num_Simpson,     // parameters for Simpson integration 1
		k_umbrella, x0,              // parameters for Simpson integration 2
		K, nyu_k, Sigma_k, pi_k,     // parameters of MG
		gammak_ij, inv_f_ik, f_i,    // responsibilities
		pi);

    Lnew=EM_L(num_sim, n, x_ij,              // # of simulation, # of snapshots, data,
	      f_i,                           // free energy differences bet. simulations
	      K, nyu_k, Sigma_k, pi_k,       // parameters of MG
	      pi);

    if (i>0) {
      dL=fabs(Lnew-L);
      ddL=dL/L*100.0;
    }
    L=Lnew;

    ++i;
    printf("# ofiteration = %3d L = %5.3e\n",i,L);
    /*    for (j=0;j<K;++j) {
      printf("nyu (%3d) = ( ",j+1);
      printf("%5.3lf ",nyu_k[j]);
      printf(") ");
      printf("Sigma (%3d) = ( ",j+1);
      printf("%5.3lf ",Sigma_k[j]);
      printf(") ");
      printf("pi=%5.3lf \n",pi_k[j]);
      }*/

    pmftempfile=efopen(pmftempfilename,"w");
    for (j=0;j<N_bin;++j) {
      x=(maxx-minx)/(double)N_bin*j+minx;
      prob=0.0;
      for (k=0;k<K;++k) {
	prob+=pi_k[k]*oneD_Gaussian(x,nyu_k[k],Sigma_k[k],pi);
      }
      fprintf(pmftempfile,"%10.8lf %10.8lf %10.8lf\n",x,prob,-1.0*log(prob));
    }
    fclose(pmftempfile);

  }

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<K;++i) {
    fprintf(outputfile,"nyu_%3d  =%8.3e ",i+1,nyu_k[i]);
    fprintf(outputfile,"Sigma_%3d=%8.3e ",i+1,Sigma_k[i]);
    fprintf(outputfile,"pi_%3d   =%8.3e \n",i+1,pi_k[i]);
  }    
  fclose(outputfile);

  pmffile=efopen(pmffilename,"w");
  for (i=0;i<N_bin;++i) {
    x=(maxx-minx)/(double)N_bin*i+minx;
    prob=0.0;
    for (j=0;j<K;++j) {
      prob+=pi_k[j]*oneD_Gaussian(x,nyu_k[j],Sigma_k[j],pi);
    }
    fprintf(pmffile,"%10.8lf %10.8lf %10.8lf\n",x,prob,-1.0*log(prob));
  }
  fclose(pmffile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-minx] minimum value of integration\n");
  printf("[-maxx] maximum value of integration\n");
  printf("[-num_Simpson] number for Simpson integration\n");
  printf("[-K] number of Gaussian\n");
  printf("[-h] help \n");
  printf("%s [-h] numsim inputfilelistname outputfilename pmffilename \n",progname);
}
