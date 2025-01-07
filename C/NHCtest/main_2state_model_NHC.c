#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "NHC.h"
#include "EF.h"
#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d;
  int NVEflag=OFF;
  int numstep=100000;
  int interval=1000,intervalout=1000,intervalnc=1000,intervalflag;
  double dt=0.0001;
  double x=0.0,v,a;
  double predict[6],correct[6];
  double V;
  double K=0.0,Kv=0.0,Pv=0.0;
  double pi;

  double TB=1.0,KBT;

  double zeta,zeta_vel,zeta_acc,predict_zeta[6],correct_zeta[6];
  double Q,tau=0.1,tau2;
  double T;
  int DOF;

  double mass=1.0;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename;
  char *outputfilename;
  char *trjfilename;
  FILE *inputfile,*trjfile,*outputfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"h",0,NULL,'h'},
    {"nve",0,NULL,'e'},
    {"temp",1,NULL,'t'},
    {"dt",1,NULL,'@'},
    {"nums",1,NULL,'s'},
    {"int",1,NULL,'i'},
    {"tau",1,NULL,'?'},
    {"nchain",1,NULL,'M'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"het:@:s:i:?:M:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'e':
      NVEflag=ON;
      break;
    case 't':
      TB=atof(optarg);
      break;
    case 'o':
      tau=atof(optarg);
      break;
    case '@':
      dt=atof(optarg);
      break;
    case 'i':
      interval=atoi(optarg);
      break;
    case 's':
      numstep=atoi(optarg);
      break;
    case '?':
      tau=atof(optarg);
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
  outputfilename = *argv;
  trjfilename    = *++argv;

  pi=acos(-1.0);

  KBT=TB;

  for (i=0;i<6;++i) correct[i] = 0.0;
  correct[0]=x;

  tau=tau/2.0/pi;
  tau2=tau*tau;
  Q=(tau2)*KBT;
  
  GearsConstant[0] = 3.0/16.0;  GearsConstant[1] = 251.0/360.0;  GearsConstant[2] = 1.0;
  GearsConstant[3] = 11.0/18.0;  GearsConstant[4] = 1.0/6.0;  GearsConstant[5] = 1.0/60.0;
  
  for (i=0;i<6;++i){
    for (j=0;j<6;++j){
      if (i != j) Telar_Matrix[i][j] = 0.0;
      else Telar_Matrix[i][i] = 1.0;
    }
  }
  
  Telar_Matrix[0][1] = 1.0;
  Telar_Matrix[0][2] = 1.0;
  Telar_Matrix[0][3] = 1.0;
  Telar_Matrix[0][4] = 1.0;
  Telar_Matrix[0][5] = 1.0;
  Telar_Matrix[1][2] = 2.0;
  Telar_Matrix[1][3] = 3.0;
  Telar_Matrix[1][4] = 4.0;
  Telar_Matrix[1][5] = 5.0;
  Telar_Matrix[2][3] = 3.0;
  Telar_Matrix[2][4] = 6.0;
  Telar_Matrix[2][5] = 10.0;
  Telar_Matrix[3][4] = 4.0;
  Telar_Matrix[3][5] = 10.0;
  Telar_Matrix[4][5] = 5.0;
  
  for (j=0;j<6;++j) correct_zeta[j]=0.0;

  outputfile=efopen(outputfilename,"w");
  trjfile=efopen(trjfilename,"w");

  for (i=0;i<numstep;++i) {
    for (j=0;j<6;++j) {
      predict[j] = 0.0;
      predict_zeta[j] = 0.0;
    }
    for (j=0;j<6;++j) {
      for (k=0;k<6;++k) { 
	predict[j] += Telar_Matrix[j][k]*correct[k];
	predict_zeta[j] += Telar_Matrix[j][k]*correct_zeta[k];
      }
    }

    x=predict[0];
    v=predict[1]/dt;

    if ( x > 20.0) x-=20.0;
    else if ( x < -20.0 ) x+=20.0;      

    if ( x < -1.0 || x > 1.0 )
      V=0.0;
    else
      V=0.5*x*x-0.5;
    K=0.5*v*v;
    T=K/2.0;

    a=-x;
    if (NVEflag==OFF) {
      a-=v*zeta_vel;
      zeta_acc = 1.0/(tau2)*(T/TB-1.0);
    }

    for (j=0;j<6;++j) {
      correct[j] = predict[j]+GearsConstant[j]*(0.5*dt*dt*a-predict[2]);
      correct_zeta[j] = predict_zeta[j]+GearsConstant[j]*(0.5*dt*dt*zeta_acc-predict_zeta[2]);
    }

    x=correct[0];
    v=correct[1]/dt;
    if ( x > 20.0) x-=20.0;
    else if ( x < -20.0 ) x+=20.0;      

    zeta=correct_zeta[0];
    zeta_vel=correct_zeta[1]/dt;

    if ( x < -1.0 || x > 1.0 )
      V=0.0;
    else
      V=0.5*x*x-0.5;
    K=0.5*v*v;
    T=K/2.0;

    Kv = 0.5*Q*zeta_vel*zeta_vel;
    Pv = KBT*zeta;

    if (i%interval==0) {
      fprintf(outputfile,"%d %e %e %e %e %e %e %e\n",i,V,K,Kv,Pv,V+K+Pv+Kv,zeta,T);
      fprintf(trjfile,"%e \n",x);
    }
  }
  fclose(outputfile);
  fclose(trjfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] outputfilename trjfilename\n",progname);
}


