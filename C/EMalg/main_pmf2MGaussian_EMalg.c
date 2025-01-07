
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <netcdf.h>
#include <getopt.h>

#include "EMalg.h"
#include "K_means.h"
#include "Gaussian.h"
#include "EF.h"

#define ON 0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,n;
  double fx,fx_old,fy_old,fy,fp;

  int N,K;
  int Nx,Ny;

  int dim;
  int flag,floatflag,flagnum,flagpmf=ON;

  //int **hist;
  double **prob,*prob2,sum;

  double *x;
  double *ax;

  int *num_k;
  double **nyu_k,***Sigma_k, *pi_k,***gamma_nk, **ave_k;
  double ***dist,J=0.0,J_new,dJ;
  double lnrou,lnrou_new,dlnrou,threhold=0.000001;

  double pi,x_map[2],minx, maxx, dx=0.01, miny, maxy, dy=0.01,minpmf;

  char *inputfilename,*outputfilename,*logfilename;
  FILE *inputfile,*outputfile,*logfile/*histfile*/;
  
  char *line;
  size_t len=0;
  
  int c,d;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;
  
  pi=acos(-1.0);
  
  struct option long_opt[] = {
    {"pb",0,NULL,'p'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"hp",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'p':
      flagpmf=OFF;
      break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);   exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < 3/*4*/) {
    USAGE(progname);
    exit(1);
  }
  K = atoi(*argv);
  //  N = atoi(*++argv);
  inputfilename = *++argv;
  outputfilename = *++argv;

  inputfile=efopen(inputfilename,"r");

  prob2=(double *)gcemalloc(sizeof(double)*1);
  i=0;
  Nx=1;
  Ny=0;
  flag=OFF;
  while ((c=fscanf(inputfile,"%lf ",&fx)) != EOF) {
    if (i==0) minx=fx;
    
    fscanf(inputfile,"%lf ",&fy);

    if (fx != fx_old && i!=0) {
      if (flag==OFF)
	maxy=fy_old;
      ++Nx;
      flag=ON;
    }

    if (flag==OFF) {
      if (i==0)
	miny=fy;
      ++Ny;
    }

    fp=0.0;
    floatflag=OFF;
    while ((c=fgetc(inputfile))!='\n') {
      if (c >= '0' && c <= '9') {
	flagnum=ON;
	if (floatflag==OFF)
	  fp=fp*10.0+c-'0';
	else  {
	  dim+=1;
	  fp+=(c-'0')/pow(10.0,dim);
	}
      }
      else if (c=='.') {
	floatflag=ON;
	dim=1;
      }
      else {
	flagnum=OFF;
      }
    }     

    if (flagnum==ON)
      if (flagpmf==ON)
	prob2[i]=exp(-1.0*fp);
      else
	prob2[i]=fp;
    else 
      prob2[i]=0.0;

    fx_old=fx;
    fy_old=fy;
    ++i;
    prob2=(double *)gcerealloc(prob2,sizeof(double)*i);
  }
  maxx=fx;

  dx=(maxx-minx)/(Nx-1);
  dy=(maxy-miny)/(Ny-1);

  /**********************/
  /* minx=0;	        */
  /* miny=0;	        */
  /* maxx=10.0;	        */
  /* maxy=10.0;	        */
  /* 		        */
  /* dx=(maxx-minx)/Nx; */
  /* dy=(maxy-miny)/Ny; */
  /**********************/

  fclose(inputfile);

  //hist=(int **)gcemalloc(sizeof(int *)*Nx);
  prob=(double **)gcemalloc(sizeof(double *)*Nx);
  for (i=0;i<Nx;++i) {
    //  hist[i]=(int *)gcemalloc(sizeof(int)*Ny);
    prob[i]=(double *)gcemalloc(sizeof(double)*Ny);
  }

  dist=(double ***)gcemalloc(sizeof(double **)*Nx);
  for (i=0;i<Nx;++i) {
    dist[i]=(double **)gcemalloc(sizeof(double *)*Nx);
    for (j=0;j<Ny;++j) dist[i][j]=(double *)gcemalloc(sizeof(double)*Ny);
  }

  n=0;
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {
      prob[i][j]=prob2[n];
      ++n;
    }
  }

  sum=0.0;
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {
      sum+=prob[i][j];
    }
  }
  for (i=0;i<Nx;++i) {
    for (j=0;j<Ny;++j) {
      prob[i][j]=prob[i][j]/sum;
    }
  }

  /*************************************************************************************************/
  /* for (i=0;i<Nx;++i) {									   */
  /*   for (j=0;j<Ny;++j) {									   */
  /*     hist[i][j]=(int)(N*prob[i][j]);							   */
  /*   }											   */
  /* }												   */
  /* 												   */
  /* histfile=efopen("hist.txt","w");								   */
  /* for (i=0;i<Nx;++i) {									   */
  /*   for (j=0;j<Ny;++j) {									   */
  /*     for (k=0;k<hist[i][j];++k) {								   */
  /* 	//      fprintf(histfile,"%8.3lf %8.3lf %10d\n",minx+dx*i,miny+dy*j,hist[i][j]);	   */
  /* 	fprintf(histfile,"%8.3lf %8.3lf\n",minx+dx*i,miny+dy*j);				   */
  /*     }											   */
  /*   }											   */
  /* }												   */
  /* fclose(histfile);										   */
  /*************************************************************************************************/

  /**************************/
  /* N=0;		    */
  /* for (i=0;i<Nx;++i) {   */
  /*   for (j=0;j<Ny;++j) { */
  /*     N+=hist[i][j];	    */
  /*   }		    */
  /* }			    */
  /**************************/

  nyu_k=(double **)gcemalloc(sizeof(double *)*K);
  for (i=0;i<K;++i) nyu_k[i]=(double *)gcemalloc(sizeof(double)*2);

  Sigma_k=(double ***)gcemalloc(sizeof(double **)*K);
  for (i=0;i<K;++i) {
    Sigma_k[i]=(double **)gcemalloc(sizeof(double *)*2);
    for (j=0;j<2;++j) Sigma_k[i][j]=(double *)gcemalloc(sizeof(double )*2);
  }

  ave_k=(double **)gcemalloc(sizeof(double *)*K);
  for (i=0;i<K;++i) {
    ave_k[i]=(double *)gcemalloc(sizeof(double)*2);
  }

  pi_k=(double *)gcemalloc(sizeof(double )*K);

  gamma_nk=(double ***)gcemalloc(sizeof(double **)*Nx);
  for (i=0;i<Nx;++i) gamma_nk[i]=(double **)gcemalloc(sizeof(double *)*Ny);
  for (i=0;i<Nx;++i) 
    for (j=0;j<Ny;++j)
      gamma_nk[i][j]=(double *)gcemalloc(sizeof(double)*K);

  /*********************************************/
  /* ax=(double *)gcemalloc(sizeof(double)*2); */
  /* ax[0]=0.0;				       */
  /* ax[1]=0.0;				       */
  /* for (i=0;i<Nx;++i) {		       */
  /*   for (j=0;j<Ny;++j) {		       */
  /*     ax[0]+=prob[i][j]*(minx+dx*i);	       */
  /*     ax[1]+=prob[i][j]*(miny+dy*i);	       */
  /*   }				       */
  /* }					       */
  /*********************************************/

  for (i=0;i<K;++i) {
    nyu_k[i][0]=/*ax[0]*/minx+dx*i;
    nyu_k[i][1]=/*ax[1]*/miny+dy*i;
  }

  x=(double *)gcemalloc(sizeof(double)*2);

  for (i=0;i<Nx;++i) {
    x[0]=minx+i*dx;
    for (j=0;j<Ny;++j) {
      x[1]=miny+i*dy;
      for (k=0;k<K;++k) {
	dist[i][j][k]=(x[0]-nyu_k[k][0])*(x[0]-nyu_k[k][0])+(x[1]-nyu_k[k][1])*(x[1]-nyu_k[k][1]);
      }
    }
  }
  Kmeans_Estep_fprob(Nx,minx,dx,Ny,miny,dy,
		     K,prob,nyu_k,gamma_nk,dist);
  J=Kmeans_J_fprob(Nx,minx,dx,Ny,miny,dy,
	     K,prob,nyu_k,gamma_nk,dist);

  i=0;

  printf("# ofiteration = %3d J = %8.3lf\n",i,J);
  for (j=0;j<K;++j) {
    printf("nyu (%3d) = ( ",j);
    for (k=0;k<2;++k) {
      printf("%8.3lf ",nyu_k[j][k]);
    }
    printf(")\n");
  }

  while (dJ > threhold || i==0) {

    Kmeans_Estep_fprob(Nx,minx,dx,Ny,miny,dy,
		       K,prob,nyu_k,gamma_nk,dist);

    Kmeans_Mstep_fprob(Nx,minx,dx,Ny,miny,dy,
		       K,prob,nyu_k,gamma_nk);
    
    J_new=Kmeans_J_fprob(Nx,minx,dx,Ny,miny,dy,
			 K,prob,nyu_k,gamma_nk,dist);

    dJ=fabs(J_new-J);
    J=J_new;

    ++i;
    printf("# ofiteration = %3d J = %8.3lf\n",i,J);
    for (j=0;j<K;++j) {
      printf("nyu (%3d) = ( ",j);
      for (k=0;k<2;++k) {
	printf("%8.3lf ",nyu_k[j][k]);
      }
      printf(")\n");
    }
  }

  num_k=(int *)gcemalloc(sizeof(int)*K);

  for (i=0;i<K;++i) {
    Sigma_k[i][0][0]=0.0;
    Sigma_k[i][1][1]=0.0;    
    ave_k[i][0]=0.0;
    ave_k[i][1]=0.0;    
  }

  for (i=0;i<Nx;++i) {
    x[0]=minx+i*dx;
    for (j=0;j<Ny;++j) {
      x[1]=miny+j*dy;
      for (k=0;k<K;++k) {
	//	pi_k[k]=gamma_nk[i][j][k];
	Sigma_k[k][0][0]+=gamma_nk[i][j][k]*x[0]*x[0]*prob[i][j];
	Sigma_k[k][1][1]+=gamma_nk[i][j][k]*x[1]*x[1]*prob[i][j];    
	ave_k[k][0]+=gamma_nk[i][j][k]*x[0]*prob[i][j];
	ave_k[k][1]+=gamma_nk[i][j][k]*x[1]*prob[i][j];    
      }
    }
  }

  for (i=0;i<K;++i) {
    pi_k[i]=1.0;
  }

  for (i=0;i<K;++i) {
    Sigma_k[i][0][0]-=ave_k[i][0]*ave_k[i][0];
    Sigma_k[i][1][1]-=ave_k[i][1]*ave_k[i][1];    
  }

  for (i=0;i<K;++i) {
    Sigma_k[i][0][0]=/*2000.0*/1000.0/*100.0*//*1.0*/;
    Sigma_k[i][1][0]=0.0;
    Sigma_k[i][0][1]=0.0;
    Sigma_k[i][1][1]=/*2000.0*/1000.0/*100.0*//*1.0*/;
  }

  lnrou=EM_lnrou_fprob(Nx,minx,dx,Ny,miny,dy,
		       K,prob,nyu_k,Sigma_k,pi_k,pi);
  i=0;
  dlnrou=0.0;
  printf("# ofiteration = %3d ln(rou) = %8.3lf\n",i,lnrou);
  for (j=0;j<K;++j) {
    printf("nyu (%3d) = ( ",j);
    for (k=0;k<2;++k) {
      printf("%8.3lf ",nyu_k[j][k]);
    }
    printf(")\n");
  }

  while (dlnrou > threhold || i==0) {

    E_step_fprob(Nx,minx,dx,Ny,miny,dy,K,prob,nyu_k,Sigma_k,pi_k,gamma_nk,pi);

    M_step_fprob(Nx,minx,dx,Ny,miny,dy,K,prob,nyu_k,Sigma_k,pi_k,gamma_nk);
    
    lnrou_new=EM_lnrou_fprob(Nx,minx,dx,Ny,miny,dy,K,prob,nyu_k,Sigma_k,pi_k,pi);

    dlnrou=fabs(lnrou_new-lnrou);
    lnrou=lnrou_new;

    ++i;
    printf("# ofiteration = %3d ln(rou) = %8.3lf\n",i,lnrou);
    for (j=0;j<K;++j) {
      printf("nyu (%3d) = ( ",j);
      for (k=0;k<2;++k) {
	printf("%8.3lf ",nyu_k[j][k]);
      }
      printf(")\n");
    }
  }

  prob2=Create_mixed_twoD_GaussianMap(minx, maxx, dx,
				      miny, maxy, dy,
				      nyu_k,  Sigma_k, pi_k, K,  pi);

  if (flagpmf==ON) {

    n=0;
    for (i=0;i<Nx;++i) {
      for (j=0;j<Ny;++j){
	prob2[n]=-1.0*log(prob2[n]);
	++n;
      }
    }

    n=0;
    minpmf=prob2[0];
    for (i=0;i<Nx;++i) {
      for (j=0;j<Ny;++j){
	if (minpmf>prob2[n])
	  minpmf=prob2[n];
	++n;
      }
    }
    
    n=0;
    for (i=0;i<Nx;++i) {
      for (j=0;j<Ny;++j){
	prob2[n]=prob2[n]-minpmf;
	++n;
      }
    }
  }

  outputfile=efopen(outputfilename,"w");
  n=0;
  x_map[0]=minx;
  for (i=0;x_map[0]<maxx;++i) {
    x_map[0]=minx+dx*i;
    x_map[1]=miny;
    for (j=0;x_map[1]<maxy;++j){
      x_map[1]=miny+dy*j;	
      fprintf(outputfile,"%8.3lf %8.3lf %8.3e\n",x_map[0],x_map[1],prob2[n]);
      ++n;
    }
  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename \n",progname);
}
