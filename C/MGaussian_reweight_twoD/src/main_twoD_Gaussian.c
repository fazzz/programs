
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <netcdf.h>
#include <getopt.h>

#include "EM_reweight_twoD.h"

#include "EM_twoD.h"

#include "Simpson_integ_twoD.h"

#include "Gaussian_twoD.h"

#include "EF.h"

#define ON 0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,num,dummy;
  double f1,f2;

  int reweight_flag=ON;

  int periodicflag=OFF;
  double periodicity;

  int num_sim,*n,n_sum;        // # of simulation, # of snapshots
  double ***x_ij;              // raw data
  double *minx, *maxx; // min and max values of data

  double **k_umbrella,**x0; // parameters for umbrella sampling
  
  // parameters of Mixed-Gaussian
  int K; // # of Gaussian 
  double **nyu_k,***Sigma_k, *pi_k;  // mean,variance,weight, # of snapshots for k-th prob.
  double ***gammak_ij, **z_ik, *z_i; // responsibilities
  double L,Lnew,dL=0.0,ddL=1.0,threhold=1.0e-4; // threhold for conversion etc.

  // parameters for Simpson integration
  int num_Simpson=10000; // # of bins (define width of integration)

  int N_bin=20;
  double x[2],prob;

  double S,sum;

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
    {"num_Simpson",1,NULL,'n'},
    {"num_bin",1,NULL,'b'},
    {"K",1,NULL,'K'},
    {"pmftemp",1,NULL,'p'},
    {"tol",1,NULL,'t'},
    {"P",1,NULL,'P'},
    {"pi",0,NULL,'i'},
    {"no_rw",0,NULL,'o'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  minx=(double *)gcemalloc(sizeof(double)*2);
  maxx=(double *)gcemalloc(sizeof(double)*2);
  
  while((c=getopt_long(argc,argv,"hiob:K:n:p:t:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'K':
      K = atoi(optarg);  break;
    case 'P':
      periodicflag = ON; periodicity=atof(optarg); break;
    case 'o':
      reweight_flag=OFF;  break;
    case 'n':
      num_Simpson=atoi(optarg);  break;
    case 'b':
      N_bin=atoi(optarg);  break;
    case 'p':
      pmftempfilename=optarg;  break;
    case 't':
      threhold=atof(optarg);  break;
    case 'i':
      periodicflag=ON; periodicity=2.0*pi; 
      minx[0]=-pi; maxx[0]=pi; minx[1]=-pi; maxx[1]=pi; break;
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

  x_ij=(double ***)gcemalloc(sizeof(double **)*2);
  for (i=0;i<2;++i) {
    x_ij[i]=(double **)gcemalloc(sizeof(double *)*num_sim);
    for (j=0;j<num_sim;++j) x_ij[i][j]=(double *)gcemalloc(sizeof(double));
  }
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
      d=fscanf(inputfile,"%lf %lf",&f1,&f2);
      x_ij[0][i]=(double *)gcerealloc(x_ij[0][i],sizeof(double)*(num+1));
      x_ij[1][i]=(double *)gcerealloc(x_ij[1][i],sizeof(double)*(num+1));
      if (periodicflag==ON) {
	while (f1<minx[0]) {
	  f1+=periodicity;
	}
	while (f1>maxx[0]) {
	  f1-=periodicity;
	}
      }
      if (periodicflag==ON) {
	while (f2<minx[1]) {
	  f2+=periodicity;
	}
	while (f2>maxx[1]) {
	  f2-=periodicity;
	}
      }
      x_ij[0][i][num]=f1;
      x_ij[1][i][num]=f2;
      ++num;
    }
    fclose(inputfile);
    n[i]=num-1;
    n_sum+=n[i];
  }
  fclose(inputfilelist);

  if (periodicflag==OFF) {
    for (i=0;i<2;++i) {
      minx[i] = x_ij[i][0][0];
      maxx[i] = x_ij[i][0][0];
      for (j=0;j<num_sim;++j) {
	for (k=0;k<n[j];++k) {
	  if (minx[i] > x_ij[i][j][k]) minx[i]=x_ij[i][j][k];
	  if (maxx[i] < x_ij[i][j][k]) maxx[i]=x_ij[i][j][k];
	}
      }
    }
  }

  k_umbrella=(double **)gcemalloc(sizeof(double *)*num_sim);
  for (i=0;i<num_sim;++i) k_umbrella[i]=(double *)gcemalloc(sizeof(double)*2);

  x0=(double **)gcemalloc(sizeof(double *)*num_sim);
  for (i=0;i<num_sim;++i) x0[i]=(double *)gcemalloc(sizeof(double)*2);

  metadatafile=efopen(metadatafilename,"r");
  for (i=0;i<num_sim;++i) {

    fscanf(metadatafile,"%lf",&x0[i][0]); 
    fscanf(metadatafile,"%lf",&k_umbrella[i][0]); 

    fscanf(metadatafile,"%lf",&x0[i][1]); 
    fscanf(metadatafile,"%lf",&k_umbrella[i][1]); 

    if (periodicflag==ON) {
      while (x0[i][0]<minx[0]) {
	x0[i][0]+=periodicity;
      }
      while (x0[i][0]>maxx[0]) {
	x0[i][0]-=periodicity;
      }
    }
    if (periodicflag==ON) {
      while (x0[i][1]<minx[1]) {
	x0[i][1]+=periodicity;
      }
      while (x0[i][1]>maxx[1]) {
	x0[i][1]-=periodicity;
      }
    }

  }
  fclose(metadatafile);
  
  nyu_k=(double **)gcemalloc(sizeof(double *)*K);
  for (i=0;i<K;++i) nyu_k[i]=(double *)gcemalloc(sizeof(double)*2);

  Sigma_k=(double ***)gcemalloc(sizeof(double **)*K);
  for (i=0;i<K;++i) {
    Sigma_k[i]=(double **)gcemalloc(sizeof(double *)*2);
    for (j=0;j<2;++j) Sigma_k[i][j]=(double *)gcemalloc(sizeof(double )*2);
  }

  pi_k=(double *)gcemalloc(sizeof(double)*K);

  //Initialization
  for (i=0;i<K;++i) {
    for (j=0;j<2;++j) {
      nyu_k[i][j]=((maxx[j]-minx[j])/(double)K*i+minx[j]);
    }
    Sigma_k[i][0][0]=10.0;
    Sigma_k[i][0][1]=0.0;
    Sigma_k[i][1][0]=0.0;
    Sigma_k[i][1][1]=10.0;
  }

  gammak_ij=(double ***)gcemalloc(sizeof(double **)*K);
  for (i=0;i<K;++i) {
    gammak_ij[i]=(double **)gcemalloc(sizeof(double *)*num_sim);
    for (j=0;j<num_sim;++j) {
      gammak_ij[i][j]=(double *)gcemalloc(sizeof(double)*n[j]);
    }
  }

  z_ik=(double **)gcemalloc(sizeof(double *)*num_sim);
  for (i=0;i<num_sim;++i) {
    z_ik[i]=(double *)gcemalloc(sizeof(double)*K);
  }
  z_i=(double *)gcemalloc(sizeof(double)*num_sim);

  for (k=0;k<K;++k) {
    S=Simpson_integ_twoD_Gaussian(num_Simpson,minx,maxx,nyu_k[k],Sigma_k[k],pi);
    pi_k[k]=1.0/S/K;
  }

  /***********************/
  /* for (k=0;k<K;++k) { */
  /*   pi_k[k]=1.0/K;	 */
  /* }			 */
  /***********************/

  i=0;
  dL=threhold;
  while (ddL >= threhold || i==0) {

    if (reweight_flag==ON) {
      E_step_twoD_rw(num_sim,n,x_ij,                          // # of simulation, # of snapshots, data,
    		     minx,maxx,num_Simpson,   // parameters for Simpson integration 1
    		     k_umbrella,x0,                           // parameters for Simpson integration 2
    		     K,nyu_k,Sigma_k,pi_k,                    // parameters of MG
    		     gammak_ij,z_ik,z_i,                      // responsibilities
    		     pi,periodicflag,periodicity);
    
      M_step_twoD_rw(num_sim, n, n_sum, x_ij,                 // # of simulation, # of snapshots, data,
    		     minx, maxx, num_Simpson, // parameters for Simpson integration 1
    		     k_umbrella, x0,                          // parameters for Simpson integration 2
    		     K, nyu_k, Sigma_k, pi_k,                 // parameters of MG
    		     gammak_ij, z_ik, z_i,                    // responsibilities
    		     pi,periodicflag,periodicity);
    
      Lnew=EM_L_twoD_rw(num_sim, n, x_ij,              // # of simulation, # of snapshots, data,
    			z_i,                           // free energy differences bet. simulations
    			K, nyu_k, Sigma_k, pi_k,       // parameters of MG
    			pi);
    }
    else {
      E_step_twoD(num_sim,n,x_ij,              // # of simulation, # of snapshots, data,
		  K, nyu_k, Sigma_k, pi_k,     // parameters of MG
		  gammak_ij,                   // responsibilities
		  pi);
      
      M_step_twoD(num_sim, n, n_sum, x_ij,     // # of simulation, # of snapshots, data,
		  K, nyu_k, Sigma_k, pi_k,     // parameters of MG
		  gammak_ij,                   // responsibilities
		  pi);

      Lnew=EM_L_twoD(num_sim, n, x_ij,              // # of simulation, # of snapshots, data,
		     K, nyu_k, Sigma_k, pi_k,       // parameters of MG
		     pi);
    }

    if (i>0) {
      dL=fabs(Lnew-L);
      ddL=dL/fabs(L)*100.0;
    }
    L=Lnew;

    sum=0.0;
    for (j=0;j<K;++j) {
      sum+=pi_k[j];
    }

    ++i;
    printf("%3d %5.3e\n",i,L);
    /**********************************************************************************************/
    /* for (j=0;j<K;++j) {									  */
    /*   printf("nyu_d  = %8.3e %8.3e ",j+1,nyu_k[j][0],nyu_k[j][1]);				  */
    /*   printf("Sigma_d = %8.3e %8.3e %8.3e %8.3e ",						  */
    /* 	      j+1,Sigma_k[j][0][0],Sigma_k[j][0][1],Sigma_k[j][1][0],Sigma_k[j][1][1]);		  */
    /*   printf("pi_d   = %8.3e \n",j+1,pi_k[j]);						  */
    /* }    											  */
    /**********************************************************************************************/

    pmftempfile=efopen(pmftempfilename,"w");
    for (j=0;j<N_bin;++j) {
      x[0]=(maxx[0]-minx[0])/(double)N_bin*j+minx[0];
      for (k=0;k<N_bin;++k) {
	x[1]=(maxx[1]-minx[1])/(double)N_bin*k+minx[1];
	prob=0.0;
	for (l=0;l<K;++l) {
	  prob+=pi_k[l]*twoD_Gaussian(x,nyu_k[l],Sigma_k[l],pi);
	}
	//	fprintf(pmftempfile,"%10.8lf %10.8lf %10.8lf %10.8lf\n",x[0],x[1],prob,-1.0*log(prob));
	fprintf(pmftempfile,"%10.8lf %10.8lf %10.8lf\n",x[0],x[1],-1.0*log(prob));
      }
      fprintf(pmftempfile,"\n");
    }
    fclose(pmftempfile);

  }

  outputfile=efopen(outputfilename,"w");
  for (j=0;j<K;++j) {
    fprintf(outputfile,"nyu_d  = %8.3e %8.3e ",j+1,nyu_k[j][0],nyu_k[j][1]);
    fprintf(outputfile,"Sigma_d = %8.3e %8.3e %8.3e %8.3e ",
	    j+1,Sigma_k[j][0][0],Sigma_k[j][0][1],Sigma_k[j][1][0],Sigma_k[j][1][1]);
    fprintf(outputfile,"pi_d   = %8.3e \n",j+1,pi_k[j]);
  }    
  fclose(outputfile);

  pmffile=efopen(pmffilename,"w");
  for (i=0;i<N_bin;++i) {
    x[0]=(maxx[0]-minx[0])/(double)N_bin*i+minx[0];
    for (j=0;j<N_bin;++j) {
      x[1]=(maxx[1]-minx[1])/(double)N_bin*j+minx[1];
      prob=0.0;
      for (k=0;k<K;++k) {
	prob+=pi_k[k]*twoD_Gaussian(x,nyu_k[k],Sigma_k[k],pi);
      }
      //      fprintf(pmffile,"%10.8lf %10.8lf %10.8lf %10.8lf\n",x[0],x[1],prob,-1.0*log(prob));
      fprintf(pmffile,"%10.8lf %10.8lf %10.8lf \n",x[0],x[1],-1.0*log(prob));
    }
    fprintf(pmffile,"\n");
  }
  fclose(pmffile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-num_Simpson] number for Simpson integration\n");
  printf("[-K] number of Gaussian\n");
  printf("[-h] help \n");
  printf("%s [-h] numsim inputfilelistname outputfilename pmffilename \n",progname);
}
