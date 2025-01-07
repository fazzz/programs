
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "QUA.h"
#include "EF.h"
#include "IO.h"
#include "LA.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  double u[3],x[4],q[4],x_roted[4],lenu;

  char *line;
  size_t len=0;

  char *progname;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double pi;
  double psi;

  progname=argv[0];

  pi=acos(-1.0);

  x[0]=0.0;
  /*************/
  /* x[1]=0.0; */
  /* x[2]=0.0; */
  /* x[3]=1.0; */
  /*************/

  psi=pi/2.0;
  q[0]=cos(0.5*psi);
  q[1]=0.0;
  q[2]=sin(0.5*psi);
  q[3]=0.0;

  for (i=0;i<4;++i) x_roted[i]=0.0;

  while((c=getopt(argc,argv,"hp:"))!=-1) {
    switch(c) {
    case 'h':
      USAGE(progname);
      exit(1);
    case 'p':
      psi=atof(optarg);
      psi=psi/180.0*pi;
      break;
    default:
      USAGE(progname);
      exit(1);
    }
  }

  argc-=optind;
  argv+=optind;
  
  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }
  x[1]=atof(*argv);
  x[2]=atof(*++argv);
  x[3]=atof(*++argv);
  u[0]=atof(*++argv);
  u[1]=atof(*++argv);
  u[2]=atof(*++argv);

  lenu=inprod(u,u,3);
  if (lenu==0.0) {
    printf("eror\n");
    exit(1);
  }
  for (i=0;i<3;++i) u[i]=u[i]/sqrt(lenu);

  q[0]=cos(0.5*psi);
  for (i=1;i<4;++i) q[i]=sin(0.5*psi)*u[i-1];

  qua_rot(x,q,x_roted);
  
  for (i=1;i<4;++i)
    printf("%10.4lf\n",x[i]);
  printf("\n");
  for (i=1;i<4;++i)
    printf("%10.4lf\n",x_roted[i]);
}

int USAGE(char *progname) {
  printf("USAGE:%s \n",progname);
  printf(" -h help\n");
  printf(" -p psi(degree)\n");
  printf(" x y z ux uy uz\n");
}

 
