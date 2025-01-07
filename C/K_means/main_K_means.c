
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "K_means.h"
#include "EF.h"

#define ON 1
#define OFF 0

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k;
  double f;

  int N,K;
  
  double M;
  double **x,**nyu,**gamma_nk;
  double avex[2],varx[2];
  int *vec;
  double **dist,J=0.0,J_new,dJ,threhold=0.00000001;

  char *inputfilename,*outputfilename,*logfilename;
  FILE *inputfile,*outputfile,*logfile;

  int logflag=OFF;

  char *line;
  size_t len=0;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;
  
  struct option long_opt[] = {
    {"log",1,NULL,'l'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"hl:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'l':
      logflag=ON;
      logfilename=optarg;
      break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);   exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  K = atoi(*argv);
  inputfilename = *++argv;
  outputfilename = *++argv;

  inputfile=efopen(inputfilename,"r");
  i=0;
  x=(double **)gcemalloc(sizeof(double *)*1);
  x[0]=(double *)gcemalloc(sizeof(double)*2);
  while ((c=fscanf(inputfile,"%lf",&f))!=EOF) {
    x[i][0]=f;
    if (c=fscanf(inputfile,"%lf",&f)!=EOF)
      x[i][1]=f;
    else {
      exit(1);
      printf("inputfile error ! data not match");
    }
    x=(double **)gcerealloc(x,sizeof(double *)*(i+1));
    x[i+1]=(double *)gcemalloc(sizeof(double )*2);
    ++i;
  }
  fclose(inputfile);
  N=i;

  nyu=(double **)gcemalloc(sizeof(double *)*K);
  for (i=0;i<K;++i) nyu[i]=(double *)gcemalloc(sizeof(double)*2);
  gamma_nk=(double **)gcemalloc(sizeof(double *)*N);
  for (i=0;i<N;++i) gamma_nk[i]=(double *)gcemalloc(sizeof(double)*K);

  vec=(int *)gcemalloc(sizeof(int)*N);

  dist=(double **)gcemalloc(sizeof(double *)*N);
  for (i=0;i<N;++i) dist[i]=(double *)gcemalloc(sizeof(double)*K);

  /**************************************/
  /* for (i=0;i<2;++i) {	        */
  /*   avex[i]=0.0;		        */
  /*   varx[i]=0.0;		        */
  /*   for (j=0;j<N;++j) {	        */
  /*     avex[i]+=x[j][i];	        */
  /*     varx[i]+=x[j][i]*x[j][i];      */
  /*   }			        */
  /*   avex[i]/=N;		        */
  /*   varx[i]/=N;		        */
  /*   varx[i]-=avex[i]*avex[i];        */
  /* }				        */
  /* 				        */
  /* for (i=0;i<K;++i) {	        */
  /*   if (i%2==0) {		        */
  /*     nyu[i][0]=avex[0]+varx[0]*i/2; */
  /*     nyu[i][1]=avex[1]+varx[1]*i/2; */
  /*   }			        */
  /*   else {			        */
  /*     nyu[i][0]=avex[0]-varx[0]*i/2; */
  /*     nyu[i][1]=avex[1]-varx[1]*i/2; */
  /*   }			        */
  /* }				        */
  /**************************************/

  for (i=0;i<K;++i) {
    nyu[i][0]=x[i][0];
    nyu[i][1]=x[i][1];
  }

  for (i=0;i<N;++i) {
    for (j=0;j<K;++j) {
      dist[i][j]=(x[i][0]-nyu[j][0])*(x[i][0]-nyu[j][0])
	        +(x[i][1]-nyu[j][1])*(x[i][1]-nyu[j][1]);
    }
  }
  Kmeans_Estep(N,K,x,nyu,gamma_nk,dist);
  J=Kmeans_J(N,K,x,nyu,gamma_nk,dist);

  i=0;
  dJ=0.0;
  printf("# ofiteration = %3d J = %8.3lf\n",i,J);
  for (j=0;j<K;++j) {
    printf("nyu (%3d) = ( ",j);
    for (k=0;k<2;++k) {
      printf("%8.3lf ",nyu[j][k]);
    }
    printf(")\n");
  }

  while (dJ > threhold || i==0) {

    Kmeans_Estep(N,K,x,nyu,gamma_nk,dist);

    Kmeans_Mstep(N,K,x,nyu,gamma_nk);
    
    J_new=Kmeans_J(N,K,x,nyu,gamma_nk,dist);

    dJ=fabs(J_new-J);
    J=J_new;

    ++i;
    printf("# ofiteration = %3d J = %8.3lf\n",i,J);
    for (j=0;j<K;++j) {
      printf("nyu (%3d) = ( ",j);
      for (k=0;k<2;++k) {
	printf("%8.3lf ",nyu[j][k]);
      }
      printf(")\n");
    }
  }

  /***************************************************/
  /* printf("# ofiteration = %3d J = %8.3lf\n",i,J); */
  /* 						     */
  /* for (i=0;i<K;++i) {			     */
  /*   printf("nyu (%3d) = ( ",i);		     */
  /*   for (j=0;j<2;++j) {			     */
  /*     printf("%8.3lf ",nyu[i][j]);		     */
  /*   }					     */
  /*   printf(")\n");				     */
  /* }						     */
  /***************************************************/

  for (i=0;i<N;++i) {
    for (j=0;j<K;++j) {
      if (gamma_nk[i][j]!=0.0)
	vec[i]=j;
    }
  }

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<N;++i) {
    fprintf(outputfile,"%8.3lf %8.3lf %3d\n",x[i][0],x[i][1],vec[i]);
  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename \n",progname);
}


