
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <math.h>
#include <ctype.h>
#include <netcdf.h>
#include <getopt.h>

#include "Optimize_BFGS_GBRW_multiD_2.h"
#include "Simpson_integ_GBRW_multiD_2.h"
#include "EF.h"

#define ON  0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,ii,jj,num,numnum,dummy;
  int dim=2;
  int K_2;
  double f,dx,dy,lnS;

  double k_B=1.98723e-3;             
  double UNITT=418.4070;
  double T=300;

  int interval=1;

  struct data dat;

  int n_sum;   
  double *minx, *maxx;

  //  double **w_l_k,**g_l_k;
  double *w_k,*g_k;
  
  int N_bin=/*400*/40;
  double *x,prob,p;

  double S,sum;

  int periodicflag=ON;
  int dismissflag=OFF;
  double delta,delta2;
  double periodicity,*minv,*maxv;

  double **bk_l_i;

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
    {"dim",1,NULL,'d'},
    {"interval",1,NULL,'m'},
    {"bin",1,NULL,'b'},
    {"num_Simpson",1,NULL,'n'},
    {"dissmiss",0,NULL,'s'},
    {"T",1,NULL,'t'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"hsd:n:m:t:b:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'd':
      dim=atoi(optarg); break;
    case 'b':
      N_bin=atoi(optarg); break;
    case 'n':
      num_Simpson=atoi(optarg);  break;
    case 't':
      T=atof(optarg);     break;
    case 'm':
      interval=atoi(optarg);     break;
    case 's':
      dismissflag=ON;  break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);  exit(1);
    }
  }

  progname=*argv;  argc-=optind;  argv+=optind;
  
  dat.dim=dim;

  if (argc < 7) {
    USAGE(progname);
    exit(1);
  }
  dat.k_each = atoi(*argv);
  dat.K=1.0;
  for (i=0;i<dat.dim;++i) dat.K = dat.K*dat.k_each;
  dat.h_l = (double *)gcemalloc(sizeof(double)*dat.dim);
  for (i=0;i<dat.dim;++i) dat.h_l[i] = atof(*++argv);
  dat.n_sim         = atoi(*++argv);
  inputfilelistname = *++argv;
  metadatafilename  = *++argv;
  outputfilename    = *++argv;
  pmffilename       = *++argv;

  periodicflag=ON; periodicity=2.0*pi;
  minv=(double *)gcemalloc(sizeof(double)*dat.dim);
  maxv=(double *)gcemalloc(sizeof(double)*dat.dim);
  for (i=0;i<dat.dim;++i) {
    minv[i]=-1.0*pi; maxv[i]=1.0*pi;
  }

  dat.beta=1.0/(k_B*T/**UNITT*/);

  dat.x_l_ij=(double ***)gcemalloc(sizeof(double **)*dat.dim);
  for (i=0;i<dat.dim;++i) dat.x_l_ij[i]=(double **)gcemalloc(sizeof(double *)*(dat.n_sim));
  for (i=0;i<dat.dim;++i) for (j=0;j<dat.n_sim;++j) dat.x_l_ij[i][j]=(double *)gcemalloc(sizeof(double));

  dat.n=(int *)gcemalloc(sizeof(int)*(dat.n_sim));

  n_sum=0;
  inputfilelist=efopen(inputfilelistname,"r");
  for (i=0;i<dat.n_sim;++i) {
    getline(&line,&len,inputfilelist);
    l=strlen(line);
    line[l-1]='\0';

    inputfile=efopen(line,"r");
    if ( dismissflag==ON ) getline(&line,&len,inputfile);

    num=0; numnum=0;
    d = 1;
    while ( d != -1  )  {
      fscanf(inputfile,"%d",&dummy);
      for (j=0;j<dat.dim;++j) {
	d=fscanf(inputfile,"%lf",&f);
	if (numnum%interval==0) {
	  dat.x_l_ij[j][i]=(double *)gcerealloc(dat.x_l_ij[j][i],sizeof(double)*(num+1));
	  if (periodicflag==ON) {
	    while (f<minv[j]) {
	      f+=periodicity;
	    }
	    while (f>maxv[j]) {
	      f-=periodicity;
	    }
	  }
	  dat.x_l_ij[j][i][num]=f;
	  if (j==dat.dim-1)  ++num;
	}
      }
      ++numnum;
    }
    fclose(inputfile);
    dat.n[i]=num-1;
    n_sum+=dat.n[i];
  }
  fclose(inputfilelist);

  dat.k_l_i=(double **)gcemalloc(sizeof(double *)*dat.dim);
  for (i=0;i<dat.dim;++i) dat.k_l_i[i]=(double *)gcemalloc(sizeof(double)*dat.n_sim);

  dat.x_l_i=(double **)gcemalloc(sizeof(double *)*dat.dim);
  for (i=0;i<dat.dim;++i) dat.x_l_i[i]=(double *)gcemalloc(sizeof(double)*dat.n_sim);

  metadatafile=efopen(metadatafilename,"r");
  for (i=0;i<dat.n_sim;++i) {
    for (j=0;j<dat.dim;++j) {
      fscanf(metadatafile,"%lf",&(dat.x_l_i[j][i]));
      if (periodicflag==ON) {
	while (dat.x_l_i[j][i]<minv[j]) {
	  dat.x_l_i[j][i]+=periodicity;
	}
	while (dat.x_l_i[j][i]>maxv[j]) {
	  dat.x_l_i[j][i]-=periodicity;
	}
      }
    }
    for (j=0;j<dat.dim;++j) {
      fscanf(metadatafile,"%lf",&(dat.k_l_i[j][i])); 
    }
  }
  fclose(metadatafile);
  
  minx = (double *)gcemalloc(sizeof(double)*dat.dim);
  maxx = (double *)gcemalloc(sizeof(double)*dat.dim);
  for (k=0;k<dat.dim;++k) {
    minx[k] = dat.x_l_ij[k][0][0];
    maxx[k] = dat.x_l_ij[k][0][0];
    for (i=0;i<dat.n_sim;++i) {
      for (j=0;j<dat.n[i];++j) {
	if (minx[k] > dat.x_l_ij[k][i][j]) minx[k]=dat.x_l_ij[k][i][j];
	if (maxx[k] < dat.x_l_ij[k][i][j]) maxx[k]=dat.x_l_ij[k][i][j];
      }
    }
  }

  dat.x_l_k=(double **)gcemalloc(sizeof(double *)*dat.dim);
  for (i=0;i<dat.dim;++i) {
    dat.x_l_k[i]=(double *)gcemalloc(sizeof(double)*dat.K);
  }
  if (dat.dim==1) {
    dx=(maxx[0]-minx[0])/(dat.k_each);
    for (i=0;i<dat.k_each;++i) {
      dat.x_l_k[0][i]=minx[0]+dx*i;
    }
  }
  if (dat.dim==2) {
    dx=(maxx[0]-minx[j])/(dat.k_each);
    dy=(maxx[j]-minx[j])/(dat.k_each);
    k=0;
    for (i=0;i<dat.k_each;++i) {
      for (j=0;j<dat.k_each;++j) {
	  dat.x_l_k[0][k]=minx[0]+dx*i;
	  dat.x_l_k[1][k]=minx[1]+dy*j;
	  ++k;
      }
    }
  }
  
  /********************************************************************************/
  /* w_l_k=(double **)gcemalloc(sizeof(double *)*dat.dim);			  */
  /* for (i=0;i<dat.dim;++i)  w_l_k[i]=(double *)gcemalloc(sizeof(double)*dat.K); */
  /* g_l_k=(double **)gcemalloc(sizeof(double *)*dat.dim);			  */
  /* for (i=0;i<dat.dim;++i)  g_l_k[i]=(double *)gcemalloc(sizeof(double)*dat.K); */
  /********************************************************************************/

  w_k=(double *)gcemalloc(sizeof(double)*dat.K);
  g_k=(double *)gcemalloc(sizeof(double)*dat.K);

  dat.C_l_ik=(double ***)gcemalloc(sizeof(double **)*dat.dim);
  for (i=0;i<dat.dim;++i) dat.C_l_ik[i]=(double **)gcemalloc(sizeof(double *)*dat.n_sim);
  for (i=0;i<dat.dim;++i) for (j=0;j<dat.n_sim;++j) dat.C_l_ik[i][j]=(double *)gcemalloc(sizeof(double)*dat.K);

  dat.S_l_k_ij=(double ****)gcemalloc(sizeof(double ***)*dat.dim);
  for (i=0;i<dat.dim;++i) 
    dat.S_l_k_ij[i]=(double ***)gcemalloc(sizeof(double **)*dat.K);
  for (i=0;i<dat.dim;++i) 
    for (j=0;j<dat.K;++j) 
      dat.S_l_k_ij[i][j]=(double **)gcemalloc(sizeof(double *)*dat.n_sim);
  for (i=0;i<dat.dim;++i) 
    for (j=0;j<dat.K;++j) 
      for (k=0;k<dat.n_sim;++k) 
	dat.S_l_k_ij[i][j][k]=(double *)gcemalloc(sizeof(double)*dat.n[j]);

  bk_l_i=(double **)gcemalloc(sizeof(double *)*dat.dim);
  for (l=0;l<dat.dim;++l) bk_l_i[l]=(double *)gcemalloc(sizeof(double)*dat.n_sim);
  for (l=0;l<dat.dim;++l) for (i=0;i<dat.n_sim;++i) bk_l_i[l][i]=dat.beta*dat.k_l_i[l][i];

  for (l=0;l<dim;++l) {
    for (i=0;i<dat.n_sim;++i) {
      for (k=0;k<dat.K;++k) {
	dat.C_l_ik[l][i][k]=Simpson_integ_oneD_GBRW_C_GB(num_Simpson,
							 minv[l],maxv[l],
							 dat.x_l_i[l][i],dat.x_l_k[l][k],
							 bk_l_i[l][i],dat.h_l[l]);
      }
    }
  }

  dat.A_l_k=(double **)gcemalloc(sizeof(double)*dat.dim);
  for (i=0;i<dat.dim;++i) dat.A_l_k[i]=(double *)gcemalloc(sizeof(double)*dat.K);

  for (l=0;l<dim;++l) {
    for (k=0;k<dat.K;++k) {
      dat.A_l_k[l][k]=Simpson_integ_oneD_GBRW_GB(num_Simpson,minv[l],maxv[l],dat.x_l_k[l][k],dat.h_l[l]);
      dat.A_l_k[l][k]=1.0/dat.A_l_k[l][k];
    }
  }

  for (l=0;l<dim;++l) {
    for (k=0;k<dat.K;++k) {
      for (i=0;i<dat.n_sim;++i) {
	for (j=0;j<dat.n[i];++j) {
	  delta=fabs(dat.x_l_ij[l][i][j]-dat.x_l_k[l][k]);
	  if (delta>pi) delta=2.0*pi-delta;

	  lnS=-0.5*delta*delta/dat.h_l[l];
	  dat.S_l_k_ij[l][k][i][j]=exp(lnS);
	}
      }
    }
  }

  /******************************************/
  /* optimize_lnL_BFGS_multiD_2(g_l_k,dat); */
  /******************************************/
  optimize_lnL_BFGS_multiD_2(g_k,dat);

  /*******************************************************************************/
  /* for (i=0;i<dat.dim;++i) for (j=0;j<dat.K;++j) w_l_k[i][j]=exp(g_l_k[i][j]); */
  /* sum=0.0; for (i=0;i<dat.dim;++i) for (j=0;j<dat.K;++j) sum+=w_l_k[i][j];	 */
  /* for (i=0;i<dat.dim;++i) for (j=0;j<dat.K;++j) w_l_k[i][j]=w_l_k[i][j]/sum;	 */
  /*******************************************************************************/

  for (i=0;i<dat.K;++i) w_k[i]=exp(g_k[i]);
  /***********************************************/
  /* sum=0.0; for (i=0;i<dat.K;++i) sum+=w_k[i]; */
  /* for (i=0;i<dat.K;++i) w_k[i]=w_k[i]/sum;	 */
  /***********************************************/

  /**************************************************************************************************************/
  /* sum=0.0; 												        */
  /* for (l=0;l<dim;++l) {										        */
  /*   for (k=0;k<dat.K;++k) {										        */
  /*     sum+=w_l_k[l][k]*dat.A_l_k[l][k]*Simpson_integ_oneD_GBRW_multiD_GB(num_Simpson,minv[l],maxv[l],        */
  /* 									 dat.x_l_k[l][k],dat.h_l[l]);	        */
  /*   }												        */
  /* }													        */
  /**************************************************************************************************************/

  outputfile=efopen(outputfilename,"w");
  /*************************************************/
  /* for (i=0;i<dat.dim;++i) {			   */
  /*   fprintf(outputfile,"  w[%2d] =  ",i+1);	   */
  /*   for (j=0;j<dat.K;++j) {			   */
  /*     fprintf(outputfile,"  %e, ",w_l_k[i][j]); */
  /*   }					   */
  /*   fprintf(outputfile," \n");		   */
  /* }						   */
  /*************************************************/

  for (i=0;i<dat.K;++i) {
    fprintf(outputfile,"  w[%2d] =  %e\n",i+1,w_k[i]);
  }
  fclose(outputfile);

  pmffile=efopen(pmffilename,"w");
  if (dat.dim==1) {
    x=(double *)gcemalloc(sizeof(double)*1);
    for (i=0;i<N_bin;++i) {
      x[0]=(maxx[0]-minx[0])/(double)N_bin*i+minx[0];
      prob=0.0;
      for (j=0;j<dat.K;++j) {
	delta=fabs(x[0]-dat.x_l_k[0][j]);
	if (fabs(delta)>pi) delta=2.0*pi-delta;
	//      p=w_k[j]*dat.A[j]*exp(-0.5*delta*delta/dat.h);
	//      prob+=p;
	//	prob+=w_l_k[0][j]*dat.A_l_k[0][j]*exp(-0.5*delta*delta/dat.h_l[0]);
	prob+=w_k[j]*dat.A_l_k[0][j]*exp(-0.5*delta*delta/dat.h_l[0]);
      }
      fprintf(pmffile,"%10.8lf %10.8lf %10.8lf\n",x[0],prob,-1.0*dat.beta*log(prob));
    }
  }
  else if (dat.dim==2) {
    x=(double *)gcemalloc(sizeof(double)*2);
    for (i=0;i<N_bin;++i) {
      x[0]=(maxx[0]-minx[0])/(double)N_bin*i+minx[0];
      for (j=0;j<N_bin;++j) {
	x[1]=(maxx[1]-minx[1])/(double)N_bin*j+minx[1];
	prob=0.0;
	for (k=0;k<dat.K;++k) {
	  delta=fabs(x[0]-dat.x_l_k[0][k]);
	  if (fabs(delta)>pi) delta=2.0*pi-delta;
	  delta2=fabs(x[1]-dat.x_l_k[1][k]);
	  if (fabs(delta2)>pi) delta2=2.0*pi-delta;
	  prob+=w_k[k]
	    *dat.A_l_k[0][k]*exp(-0.5*delta*delta/dat.h_l[0])
	    *dat.A_l_k[1][k]*exp(-0.5*delta2*delta2/dat.h_l[1]);
	}
	fprintf(pmffile,"%10.8lf %10.8lf %10.8lf\n",x[0],x[1],-1.0*dat.beta*log(prob));
      }
      fprintf(pmffile,"\n");
    }
  }
  else if (dat.dim==3) {
    ;
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
