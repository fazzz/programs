
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "EF.h"
#include "netcdf_mine.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i;
  int numstep=3.0*10e6;

  double hdt=0.025;
  double x=-0.1,p=1.0,e=0.1,f=-1.0;
  double K;

  double F,G;
  double kbT=1;
  double m=1,me=1;
  double omega=1;
  double s,ds;
  double a,b;

  double cff1,cff2,cff3;

  int interval=100;

  char *line;
  size_t len=0;

  char *progname;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *outputfilename;
  FILE *outputfile;

  //  struct my_netcdf_out_id_MCD nc_id_MCD;

  cff1=m*omega*omega;
  cff2=4.0*4.0/5.0/pow(4.0*sqrt(2.0)/5.0,4);
  cff3=pow(4.0*sqrt(2.0)/5.0,2);

  progname=argv[0];
  while((c=getopt(argc,argv,"hn:o:"))!=-1) {
    switch(c) {
    case 'n':
      numstep=atoi(optarg);
      break;
    case 'o':
      interval=atoi(optarg);
      break;
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

  if (argc < 1) {
    USAGE(progname);
    exit(1);
  }
  outputfilename=*argv;

  K=0.5*m*p*p+0.5*me*f*f;
  outputfile=efopen(outputfilename,"w");
  fprintf(outputfile,"%8.3lf%8.3lf%8.3lf%8.3lf%8.3lf\n",x,p,e,f,K);
  for (i=0;i<numstep;++i) {
    F=-cff1*x;
    G=-cff2*(e*e-cff3)*e;

    a=(F*p/m+G*f/me)/(2.0*kbT);
    b=(F*F/m+G*G/me)/(2.0*kbT);

    s=a/b*(cosh(hdt*sqrt(b))-1)+1.0/sqrt(b)*sinh(hdt*sqrt(b));
    ds=a/sqrt(b)*sinh(hdt*sqrt(b))+cosh(hdt*sqrt(b));

    p=(p+F*s)/ds;
    f=(f+G*s)/ds;

    x+=2.0*hdt*p;
    e+=2.0*hdt*f;

    F=-cff1*x;
    G=-cff2*(e*e-cff3)*e;

    a=(F*p/m+G*f/me)/(2.0*kbT);
    b=(F*F/m+G*G/me)/(2.0*kbT);

    s=a/b*(cosh(hdt*sqrt(b))-1)+1.0/sqrt(b)*sinh(hdt*sqrt(b));
    ds=a/sqrt(b)*sinh(hdt*sqrt(b))+cosh(hdt*sqrt(b));

    p=(p+F*s)/ds;
    f=(f+G*s)/ds; 

    K=0.5*m*p*p+0.5*me*f*f;

    if (i%interval==0) fprintf(outputfile,"%8.3lf%8.3lf%8.3lf%8.3lf%8.3lf\n",x,p,e,f,K);
  }
  fclose(outputfile);

}

int USAGE(char *progname) {
  printf("USAGE:%s \n",progname);
  printf("[-n numstep] numstep  \n");
  printf("[-o interval] interval  \n");
  printf("[-h] help  \n");
  printf("%s [-n numstep] [-o interval] [-h] outputfilename \n",progname);
}

 
