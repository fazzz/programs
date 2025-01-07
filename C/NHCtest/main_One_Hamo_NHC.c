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
  double x=5.0,v,a;
  double *predict,*correct;
  double V;
  double K=0.0,Kv=0.0,Pv=0.0;
  double omega=0.01;
  double dummy;
  double pi;

  double TB=300,KBT;
  double k_B=1.98723e-3;

  double *zeta,*zeta_vel,*zeta_acc,**predict_zeta,**correct_zeta;
  double *Q_NH,tau=0.1,tau2;
  double T;
  int DOF,M=4;

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
    case 'M':
      M=atoi(optarg);
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

  predict=(double *)gcemalloc(sizeof(double)*6);
  correct=(double *)gcemalloc(sizeof(double)*6);

  zeta=(double *)gcemalloc(sizeof(double)*M);
  zeta_vel=(double *)gcemalloc(sizeof(double)*M);
  zeta_acc=(double *)gcemalloc(sizeof(double)*M);
  predict_zeta=(double **)gcemalloc(sizeof(double *)*M);
  correct_zeta=(double **)gcemalloc(sizeof(double *)*M);
  for (i=0;i<M;++i) {
    predict_zeta[i]=(double *)gcemalloc(sizeof(double)*6);
    correct_zeta[i]=(double *)gcemalloc(sizeof(double)*6);
  }
  Q_NH=(double *)gcemalloc(sizeof(double)*M);

  DOF=1;
  KBT=k_B*TB;

  for (i=0;i<6;++i) correct[i] = 0.0;
  correct[0]=x;
  NHC_set(tau,&tau2,Q_NH,M,DOF,KBT,correct_zeta);

  /****************************************************/
  /* tau=tau/2.0/pi;				      */
  /* tau2=tau*tau; 				      */
  /* Q_NH[0]=(tau2)*DOF*KBT;   			      */
  /* 						      */
  /* GearsConstant[0] = 3.0/16.0;		      */
  /* GearsConstant[1] = 251.0/360.0;		      */
  /* GearsConstant[2] = 1.0;			      */
  /* GearsConstant[3] = 11.0/18.0;		      */
  /* GearsConstant[4] = 1.0/6.0;		      */
  /* GearsConstant[5] = 1.0/60.0;		      */
  /* 						      */
  /* for (i=0;i<6;++i){				      */
  /*   for (j=0;j<6;++j){			      */
  /*     if (i != j) Telar_Matrix[i][j] = 0.0;	      */
  /*     else Telar_Matrix[i][i] = 1.0;		      */
  /*   }					      */
  /* }						      */
  /* 						      */
  /* Telar_Matrix[0][1] = 1.0;			      */
  /* Telar_Matrix[0][2] = 1.0;			      */
  /* Telar_Matrix[0][3] = 1.0;			      */
  /* Telar_Matrix[0][4] = 1.0;			      */
  /* Telar_Matrix[0][5] = 1.0;			      */
  /* Telar_Matrix[1][2] = 2.0;			      */
  /* Telar_Matrix[1][3] = 3.0;			      */
  /* Telar_Matrix[1][4] = 4.0;			      */
  /* Telar_Matrix[1][5] = 5.0;			      */
  /* Telar_Matrix[2][3] = 3.0;			      */
  /* Telar_Matrix[2][4] = 6.0;			      */
  /* Telar_Matrix[2][5] = 10.0;			      */
  /* Telar_Matrix[3][4] = 4.0;			      */
  /* Telar_Matrix[3][5] = 10.0;			      */
  /* Telar_Matrix[4][5] = 5.0;			      */
  /* 						      */
  /* for (i=0;i<M;++i) 				      */
  /*   for (j=0;j<6;++j) 			      */
  /*     correct_zeta[i][j]=0.0;		      */
  /****************************************************/
  
  outputfile=efopen(outputfilename,"w");
  trjfile=efopen(trjfilename,"w");

  for (i=0;i<numstep;++i) {
    for (j=0;j<6;++j) {
      predict[j] = 0.0;
      //      predict_zeta[0][j] = 0.0;
    }
    for (j=0;j<6;++j) {
      for (k=0;k<6;++k) { 
	predict[j] += Telar_Matrix[j][k]*correct[k];
	//	predict_zeta[0][j] += Telar_Matrix[j][k]*correct_zeta[0][k];
      }
    }

    x=predict[0];
    v=predict[1]/dt;
    /**************************************/
    /* zeta[0]=predict_zeta[0][0];	  */
    /* zeta_vel[0]=predict_zeta[0][1]/dt; */
    /**************************************/

    NHC_update_pret(zeta,zeta_vel,predict_zeta,correct_zeta,M,dt);

    V=0.5*mass*omega*omega*x*x;
    K=0.5*mass*v*v;
    T=K/(DOF*k_B)*2.0;

    a=-omega*omega*x;
    if (NVEflag==OFF) {
      a-=v*zeta_vel[0];
      NHC_solve(zeta_vel,zeta_acc,Q_NH,M,1,KBT,tau2,T,TB);
      //      zeta_acc[0] = 1.0/(tau2)*(T/TB-1.0);
    }

    for (j=0;j<6;++j) {
      correct[j] = predict[j]+GearsConstant[j]*(0.5*dt*dt*a-predict[2]);
      //      correct_zeta[0][j] = predict_zeta[0][j]+GearsConstant[j]*(0.5*dt*dt*zeta_acc[0]-predict_zeta[0][2]);
    }

    x=correct[0];
    v=correct[1]/dt;
    //    zeta[0]=correct_zeta[0][0];
    //    zeta_vel[0]=correct_zeta[0][1]/dt;
    NHC_update_cort(zeta,zeta_vel,zeta_acc,predict_zeta,correct_zeta,M,dt);

    V=0.5*mass*omega*omega*x*x;
    K=0.5*mass*v*v;
    T=K/(DOF*k_B)*2.0;
    NHC_calcKE(zeta,zeta_vel,Q_NH,M,DOF,KBT,&Pv,&Kv);
    //    Kv = 0.5*Q_NH[0]*zeta_vel[0]*zeta_vel[0];
    //    Pv = DOF*KBT*zeta[0];

    if (i%interval==0) {
      fprintf(outputfile,"%d %e %e %e %e %e %e %e\n",i,V,K,Kv,Pv,V+K+Pv+Kv,zeta[0],T);
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


