#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PTL.h"
#include "EF.h"
#include "TOPO.h"
#include "LA.h"

#define NVT 1
#define NVE 0

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  double da=0.01;
  double pi;

  double de=1.0,d=1.0,d2;

  double T0=300,T,K0,KE;
  double k_B=1.98723e-3;
  double UNITT=418.4070;
  double KT;

  double A,B,C,D;
  double p_d1_t,p_d2_t,*p_d1,*p_d2;
  double *dih_equ1,*dih_equ2;
  double e1,e2,p_MB,c1_c2,kai;

  double Kd1=1.0;
  double dihed1,dihed2;

  double x[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *parmfilename,*outputfilename,*outputfilename2;

  FILE *parmfile,*outputfile,*outputfile2;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"temp",1,NULL,'t'},
    {"da",1,NULL,'a'},
    {"de",1,NULL,'d'},
    {"dh",1,NULL,'2'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"ht:a:d:2::",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 't':
      T0=atof(optarg);
      break;
    case 'a':
      da=atof(optarg);
      break;
    case 'd':
      de=atof(optarg);
      break;
    case '2':
      d=atof(optarg);
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  outputfilename   = *argv;
  outputfilename2  = *++argv;

  pi=acos(-1.0);
  KT=k_B*T0;

  de=de*KT;
  d=d*KT;
  
  d2=d*d;

  dih_equ1=(double *)gcemalloc(sizeof(double)*2);
  dih_equ2=(double *)gcemalloc(sizeof(double)*2);

  dih_equ1[0]=-0.5*pi;
  dih_equ1[0]=0.5*pi;

  dih_equ2[1]=0.5*pi;
  dih_equ2[1]=-0.5*pi;

  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");
  for (dihed1=-pi;dihed1<pi;dihed1+=da) {
    for (dihed2=-pi;dihed2<pi;dihed2+=da) {
    
      e1=Kd1*(1.0-cos(dihed1-dih_equ1[0]))+Kd1*(1.0-cos(dihed2-dih_equ1[1]));
      e2=Kd1*(1.0-cos(dihed1-dih_equ2[0]))+Kd1*(1.0-cos(dihed2-dih_equ2[1]));

      A=0.5*(e1+e2+de);
      B=0.5*(e1-e2-de);
      C=sqrt(B*B+d2);
      D=e1-e2-de;
  
      p_MB=A-C;

      c1_c2=(e1-p_MB)/d;
      kai=log(c1_c2);
      c1_c2=d/(e2+de-p_MB);
      kai=log(c1_c2);

      c1_c2=(e1-p_MB)/d;
      kai=log((e1-p_MB)/d);
      /********************************/
      /* c1_c2=d/(e2+de-p_MB);	      */
      /* if (c1_c2<0.0) c1_c2=-c1_c2; */
      /* kai=log(c1_c2);	      */
      /********************************/

      fprintf(outputfile,"%e %e %e \n",dihed1*180/pi,dihed2*180/pi,p_MB,c1_c2,kai);
      fprintf(outputfile2,"%e %e %e \n",dihed1*180/pi,dihed2*180/pi,kai);
    }
    fprintf(outputfile," \n");
    fprintf(outputfile2," \n");
  }
  fclose(outputfile);
  fclose(outputfile2);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename refcrdfilename parmfilename outputfilename outputfilename2 trjfilename\n",progname);
}
