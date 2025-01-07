
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "LA.h"
#include "EF.h"
#include "IO.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;

  int minmn;
  int m,n;
  double *mat,*matU,*matVT,*sv;
  double f;

  char *line;
  size_t len=0;

  char *progname;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;

  progname=argv[0];

  while((c=getopt(argc,argv,"h"))!=-1) {
    switch(c) {
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

  m=2;
  n=3;
  /*******************/
  /* scanf("%d",&m); */
  /* scanf("%d",&n); */
  /*******************/
  mat=gcemalloc(sizeof(double)*m*n);
  matU=gcemalloc(sizeof(double)*m*m);
  matVT=gcemalloc(sizeof(double)*n*n);
  minmn=n;
  if (minmn > m) minmn=m;
  sv=gcemalloc(sizeof(double)*minmn);
  mat[0]=3.0;  mat[1]=1.0; mat[2]=2.0;
  mat[3]=3.0;  mat[4]=2.0; mat[5]=1.0;

  printf("P=\n");

  for (i=0;i<m;++i) {
    for (j=0;j<n;++j) {
      printf("%10.4lf ",mat[i*n+j]);
    }
    printf("\n");
  }

  printf("\n");

  /**************************************/
  /* for (i=0;i<m;++i) {	        */
  /*   for (j=0;j<n;++j) {	        */
  /*     scanf("%lf",&f);	        */
  /*     mat[i*n+j]=f;		        */
  /*     printf("%10.4lf ",mat[i*n+j]); */
  /*   }			        */
  /*   printf("\n");		        */
  /* }				        */
  /**************************************/

  svd(mat,m,n,matU,matVT,sv);


  printf("U=\n");

  for (i=0;i<m;++i) {
    for (j=0;j<m;++j) {
      printf("%10.4lf ",matU[i*m+j]);
    }
    printf("\n");
  }

  printf("\n");

  printf("VT=\n");

  for (i=0;i<n;++i) {
    for (j=0;j<n;++j) {
      printf("%10.4lf ",matVT[i*n+j]);
    }
    printf("\n");
  }

  if (n>m) m=n;

  printf("\n");

  printf("sv=\n");

  for (i=0;i<m;++i)
    printf("%10.4lf ",sv[i]);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:%s \n",progname);
  printf("-t dt    \n");
  printf("-h help  \n");
  printf("numpoint inputfilename outputfilename  \n");
}

 
