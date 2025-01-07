
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <netcdf.h>
#include <getopt.h>

#include "MoorePenlose_InverseMatrix.h"
#include "Simpson_integ_oneD_Gaussianbase.h"
#include "EF.h"

#define ON  0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,num,numnum,dummy;
  double f;

  int interval=1;

  int min_set_flag=OFF,max_set_flag=OFF;

  int periodicflag=OFF;
  double periodicity,minv,maxv;

  int num_sim,*n,n_sum;  // # of simulation, # of snapshots
  double **x_ij;         // raw data
  double minx, maxx;     // min and max values of data

  // parameters for umbrella sampling
  double *k_umbrella,*x0;

  // Bias function
  double ***C;
  
  // parameters of Mixed-Gaussian
  int K;                   // # of Gaussian 
  double *nyu_k,*Sigma_k;  // mean,variance
  double *omega_k;         // weight

  double h=0.5;

  double L,Lnew,dL=0.0,ddL=1.0,threhold=1.0e-6,lnp; // threhold for conversion etc.

  // parameters for Simpson integration
  int num_Simpson=10000;             // # of bins (define width of integration)

  // parameters for MPI
  double *Psi,*PsiT,*MPI_Psi,*MPI_PsiT;
  double *gamma_k,*gamma_prime_i,*gamma_prime_inv_i;
  double *PsiT_Psi,*Inv_Psi_PsiT;
  double *Psi_PsiT,*Inv_PsiT_Psi;

  double *f_i;

  double **Integ_C_i_GB_k;

  int N_bin=50;
  double x,prob;

  double S,sum;

  double pi;

  char *inputfilelistname,*metadatafilename,*outputfilename,*pmffilename;
  char *pmftempfilename="PMFOUT",*logfilename="log.txt";
  FILE *inputfile,*inputfilelist,*metadatafile,*outputfile,*pmffile,*pmftempfile,*logfile;
  
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
    {"pmftemp",1,NULL,'p'},
    {"log",1,NULL,'l'},
    {"interval",1,NULL,'m'},
    {"tol",1,NULL,'t'},
    {"P",1,NULL,'P'},
    {"pi",0,NULL,'i'},
    {"no_rw",0,NULL,'o'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"hib:n:p:l:m:t:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'P':
      periodicflag = ON; periodicity=atof(optarg); break;
    case 'n':
      num_Simpson=atoi(optarg);  break;
    case 'b':
      N_bin=atoi(optarg);  break;
    case 'p':
      pmftempfilename=optarg;  break;
    case 'l':
      logfilename=optarg;  break;
    case 'm':
      interval=atoi(optarg);  break;
    case 't':
      threhold=atof(optarg);  break;
    case 'i':
      periodicflag=ON; periodicity=2.0*pi; 
      minv=-1.0*pi; maxv=1.0*pi;          break;
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
    num=0; numnum=0;
    d = 1;
    while ( d != -1  )  {
      fscanf(inputfile,"%d",&dummy);
      d=fscanf(inputfile,"%lf",&f);
      if (numnum%interval==0) {
	x_ij[i]=(double *)gcerealloc(x_ij[i],sizeof(double)*(num+1));
	if (periodicflag==ON) {
	  while (f<minv) {
	    f+=periodicity;
	  }
	  while (f>maxv) {
	    f-=periodicity;
	  }
	}
	x_ij[i][num]=f;
	++num;
      }
      ++numnum;
    }
    fclose(inputfile);
    n[i]=num-1;
    n_sum+=n[i];
  }
  fclose(inputfilelist);

  k_umbrella=(double *)gcemalloc(sizeof(double)*num_sim);
  x0=(double *)gcemalloc(sizeof(double)*num_sim);
  metadatafile=efopen(metadatafilename,"r");
  for (i=0;i</*K*/num_sim;++i) {
    fscanf(metadatafile,"%lf",&x0[i]); fscanf(metadatafile,"%lf",&k_umbrella[i]); 
    if (periodicflag==ON) {
      while (x0[i]<minv) {
	x0[i]+=periodicity;
      }
      while (x0[i]>maxv) {
	x0[i]-=periodicity;
      }
    }
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

  nyu_k=(double *)gcemalloc(sizeof(double)*/*n_sum*/K);
  Sigma_k=(double *)gcemalloc(sizeof(double)*/*n_sum/*/K);
  omega_k=(double *)gcemalloc(sizeof(double)*/*n_sum*/K);

  //Initialization
  //  for (i=0;i<n_sum/*K*/;++i) {
  l=0;
  for (i=0;i<num_sim;++i) {
    for (j=0;j<n[i];++j) {
      nyu_k[/*l*/i]=/*x_ij[i][j]*/((maxx-minx)/(double)K*i+minx);
      Sigma_k[/*l*/i]=h/*(maxx-minx)/(double)K*/;
      ++l;
    }
  }

  S=0.0;
  for (k=0;k<K/*n_sum*/;++k) S+=Simpson_integ_oneD_Gaussian_base(num_Simpson,minx,maxx,nyu_k[k],Sigma_k[k]);
  for (k=0;k<K/*n_sum*/;++k) omega_k[k]=1.0/S;

  Psi=(double *)gcemalloc(sizeof(double)*n_sum*K/*n_sum*/);
  PsiT=(double *)gcemalloc(sizeof(double)*K/*n_sum*/*n_sum);
  MPI_Psi=(double *)gcemalloc(sizeof(double)*n_sum*K/*n_sum*/);
  MPI_PsiT=(double *)gcemalloc(sizeof(double)*K/*n_sum*/*n_sum);

  PsiT_Psi=(double *)gcemalloc(sizeof(double)*n_sum*n_sum);
  Inv_Psi_PsiT=(double *)gcemalloc(sizeof(double)*n_sum*n_sum);
  Psi_PsiT=(double *)gcemalloc(sizeof(double)*/*K*/n_sum*/*K*/n_sum);
  Inv_PsiT_Psi=(double *)gcemalloc(sizeof(double)*/*K*/n_sum*/*K*/n_sum);

  l=0;
  for (i=0;i<num_sim;++i) {
    for (j=0;j<n[i];++j) {
      for (k=0;k<n_sum/*K*/;++k) {
	Psi[l*n_sum/*K*/+k]=oneD_Gaussian_base(x_ij[i][j],nyu_k[k],Sigma_k[k]);
	PsiT[k*n_sum+l]=Psi[l*n_sum/*K*/+k];
      }
      ++l;
    }
  }

  /************************************************************************************/
  /* for (i=0;i<n_sum;++i) {							      */
  /*   for (j=0;j<n_sum;++j) {							      */
  /*     PsiT_Psi[i*n_sum+j]=0.0;						      */
  /*     for (k=0;k<n_sum/\*K*\/;++k) {						      */
  /* 	PsiT_Psi[i*n_sum+j]+=Psi[i*\/\*K*\/n_sum+k]*Psi[j*\/\*K*\/n_sum+k];	      */
  /*     }									      */
  /*   }									      */
  /* }										      */
  /************************************************************************************/

  /****************************************************************************/
  /* for (i=0;i</\*K*\/n_sum;++i) {					      */
  /*   for (j=0;j</\*K*\/n_sum;++j) {					      */
  /*     Psi_PsiT[i*\/\*K*\/n_sum+j]=0.0;				      */
  /*     for (k=0;k<n_sum;++k) {					      */
  /* 	Psi_PsiT[i*\/\*K*\/n_sum+j]+=Psi[i*n_sum+k]*Psi[j*n_sum+k];	      */
  /*     }								      */
  /*   }								      */
  /* }									      */
  /****************************************************************************/

  MPginvm(Psi,  MPI_Psi, n_sum, K);
  MPginvm(PsiT, MPI_PsiT, K, n_sum);
  //  mnmult(PsiT,n_sum,K,Psi,K,n_sum,PsiT_Psi);
  //  invmat(PsiT_Psi,Inv_PsiT_Psi,n_sum);
  //  mnmult(Inv_Psi_PsiT,n_sum,n_sum,PsiT,n_sum,/*K*/n_sum,MPI_Psi);

  //  mnmult(Psi,K,n_sum,PsiT,n_sum,K,Psi_PsiT);
  //  invmat(Psi_PsiT,Inv_Psi_PsiT,/*K*/n_sum);
  //  mnmult(Inv_PsiT_Psi,/*K*/n_sum,/*K*/n_sum,Psi,/*K*/n_sum,num_sim,MPI_PsiT);

  C=(double ***)gcemalloc(sizeof(double **)*num_sim);
  for (i=0;i<num_sim;++i) C[i]=(double **)gcemalloc(sizeof(double *)*num_sim);
  for (i=0;i<num_sim;++i) for (j=0;j<num_sim;++j) C[i][j]=(double *)gcemalloc(sizeof(double)*n[j]);

  for (i=0;i<num_sim;++i) 
    for (j=0;j<num_sim;++j) 
      for (k=0;k<n[j];++k)
	C[i][j][k]=C_func(k_umbrella[i],x0[i],x_ij[j][k],periodicflag,periodicity);

  Integ_C_i_GB_k=(double **)gcemalloc(sizeof(double *)*num_sim);
  for (i=0;i<num_sim;++i) Integ_C_i_GB_k[i]=(double *)gcemalloc(sizeof(double)*K);
  for (i=0;i<num_sim;++i)
    for (j=0;j</*K*/n_sum;++j)
      Integ_C_i_GB_k[i][j]=Simpson_integ_C_oneD_Gaussian_base(num_Simpson,minx,maxx,
							      k_umbrella[i],x0[i],
							      nyu_k[j],Sigma_k[j],
							      periodicflag,periodicity);

  f_i=(double *)gcemalloc(sizeof(double)*num_sim);
  for (i=0;i<num_sim;++i) f_i[i]=1.0;

  pmftempfile=efopen(pmftempfilename,"w");
  for (j=0;j<N_bin;++j) {
    x=(maxx-minx)/(double)N_bin*j+minx;
    prob=0.0;
    for (k=0;k<n_sum/*K*/;++k) {
      prob+=omega_k[k]*oneD_Gaussian_base(x,nyu_k[k],Sigma_k[k]);
    }
    fprintf(pmftempfile,"%10.8lf %10.8lf %10.8lf\n",x,prob,-1.0*log(prob));
  }
  fprintf(pmftempfile,"\n");
  fclose(pmftempfile);

  gamma_k=(double *)gcemalloc(sizeof(double)*/*K*/n_sum);
  gamma_prime_i=(double *)gcemalloc(sizeof(double)*n_sum);
  gamma_prime_inv_i=(double *)gcemalloc(sizeof(double)*n_sum);

  L=0.0;
  for (i=0;i<num_sim;++i) {
    for (j=0;j<n[i];++j) {
      lnp=0.0;
      prob=0.0;
      for (k=0;k</*K*/n_sum;++k) {
	prob+=omega_k[k]*oneD_Gaussian_base(x_ij[i][j],nyu_k[k],Sigma_k[k]);
      }
      lnp+=f_i[i]*C[i][i][j]*prob;
      lnp=log(lnp);
      L+=lnp;
    }
  }

  i=0;

  printf("%3d %5.3e\n",i,L);
  fprintf(logfile,"%3d %5.3e\n",i,L);

  /**********************************/
  /* printf("gamma_k=( ");	    */
  /* for (k=0;k<K;++k) {	    */
  /*   printf("%3.1e ",omega_k[k]); */
  /* }				    */
  /* printf(") \n");		    */
  /**********************************/

  logfile=efopen(logfilename,"w");
  dL=threhold;
  while (ddL >= threhold || i==0) {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (k=0;k</*K*/n_sum;++k) {
      gamma_k[k]=0.0; for (j=0;j<num_sim;++j) gamma_k[k]+=n[j]*f_i[j]*Integ_C_i_GB_k[j][k];
    }
    mnmult(MPI_PsiT,n_sum,/*K*/n_sum,gamma_k,/*K*/n_sum,1,gamma_prime_i);
    for (j=0;j<n_sum;++j) gamma_prime_inv_i[j]=1.0/gamma_prime_i[j];
    mnmult(MPI_Psi,/*K*/n_sum,n_sum,gamma_prime_inv_i,n_sum,1,omega_k);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    S=0.0;
    for (k=0;k</*K*/n_sum;++k) 
      S+=omega_k[k]*Simpson_integ_oneD_Gaussian_base(num_Simpson,minx,maxx,nyu_k[k],Sigma_k[k]);
    for (k=0;k</*K*/n_sum;++k) 
      omega_k[k]=1.0/S*omega_k[k];
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    for (j=0;j<num_sim;++j) {
      f_i[j]=1.0/Simpson_integ_C_P_GB(num_Simpson,minx,maxx,
				      /*K*/n_sum, omega_k,
				      k_umbrella[j],x0[j],
				      nyu_k,Sigma_k,
				      periodicflag,periodicity);
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    Lnew=0.0;
    for (j=0;j<num_sim;++j) {
      for (k=0;k<n[j];++k) {
	lnp=0.0;
	prob=0.0;
	for (l=0;l</*K*/n_sum;++l) {
	  prob+=omega_k[l]*oneD_Gaussian_base(x_ij[j][k],nyu_k[l],Sigma_k[l]);
	}
	lnp+=f_i[j]*C[j][j][k]*prob;
	lnp=log(lnp);
	Lnew+=lnp;
      }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (i>0) {
      dL=fabs(Lnew-L);
      ddL=dL/fabs(L)*100.0;
    }
    L=Lnew;

    ++i;
    printf("%3d %5.3e\n",i,L);
    fprintf(logfile,"%3d %5.3e\n",i,L);

    /**********************************/
    /* printf("gamma_k=( ");	      */
    /* for (k=0;k</\*K*\/n_sum;++k) { */
    /*   printf("%3.1e ",omega_k[k]); */
    /* }			      */
    /* printf(") \n");		      */
    /**********************************/

    sum=0.0;
    for (k=0;k</*K*/n_sum;++k) {
      S=Simpson_integ_oneD_Gaussian_base(num_Simpson,minx,maxx,nyu_k[k],Sigma_k[k]);
      sum+=omega_k[k]*S;
    }

    pmftempfile=efopen(pmftempfilename,"a");
    for (j=0;j<N_bin;++j) {
      x=(maxx-minx)/(double)N_bin*j+minx;
      prob=0.0;
      for (k=0;k</*K*/n_sum;++k) {
	prob+=omega_k[k]*oneD_Gaussian_base(x,nyu_k[k],Sigma_k[k]);
      }
      fprintf(pmftempfile,"%10.8lf %10.8lf %10.8lf\n",x,prob,-1.0*log(prob));
    }
    fprintf(pmftempfile,"\n");
    fclose(pmftempfile);

  }
  fclose(logfile);

  sum=0.0;
  for (k=0;k</*K*/n_sum;++k) {
    S=Simpson_integ_oneD_Gaussian_base(num_Simpson,minx,maxx,nyu_k[k],Sigma_k[k]);
    sum+=omega_k[k]*S;
  }

  outputfile=efopen(outputfilename,"w");
  for (i=0;i</*K*/n_sum;++i) {
    fprintf(outputfile,"nyu_%3d  =%8.3e ",i+1,nyu_k[i]);
    fprintf(outputfile,"Sigma_%3d=%8.3e ",i+1,Sigma_k[i]);
    fprintf(outputfile,"pi_%3d   =%8.3e \n",i+1,omega_k[i]);
  }    
  fclose(outputfile);

  pmffile=efopen(pmffilename,"w");
  for (i=0;i<N_bin;++i) {
    x=(maxx-minx)/(double)N_bin*i+minx;
    prob=0.0;
    for (j=0;j</*K*/n_sum;++j) {
      prob+=omega_k[j]*oneD_Gaussian_base(x,nyu_k[j],Sigma_k[j]);
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
  printf("[--pmftemp] pmf temp file name\n");
  printf("[--log] log file name\n");
  printf("[--tol] value of tolerance\n");
  printf("[-P] periodic flag\n");
  printf("[--pi] 2pi periodic flag");
  printf("[--no_rw] no reweight option\n");
  printf("[-h] help \n");
  printf("%s [-h] numsim inputfilelistname outputfilename pmffilename \n",progname);
}
