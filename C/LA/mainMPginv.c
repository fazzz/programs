
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

  int m,n;
  double *mat,*invmat;

  char *line;
  size_t len=0;

  char *progname;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

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
  n=4;
  mat=gcemalloc(sizeof(double)*m*n);
  invmat=gcemalloc(sizeof(double)*n*m);

  mat[0]=1.0;  mat[1]=2.0; mat[2]=3.0; mat[3]=4.0;
  mat[4]=2.0;  mat[5]=4.0; mat[6]=6.0; mat[7]=8.0;

  printf("A=\n");

  for (i=0;i<m;++i) {
    for (j=0;j<n;++j) {
      printf("%10.4lf ",mat[i*n+j]);
    }
    printf("\n");
  }

  printf("\n");
  printf("A+=\n");

  MPginvm(mat,invmat,m,n);


  
  for (i=0;i<n;++i) {
    for (j=0;j<m;++j) {
      printf("%10.4lf ",invmat[i*m+j]);
    }
    printf("\n");
  }

  printf("\n");

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:%s \n",progname);
  printf("-h help  \n");
  printf("--  \n");
}

 
