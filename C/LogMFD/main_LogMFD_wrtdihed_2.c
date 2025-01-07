
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EF.h"

#include "MD.h"
#include "MD_NHC_MP1996.h"

#include "LogMFD.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k;
  int N;

  int nc=1;
  double dt=0.001,dt2,wdt2[3],wdt4[3];
  double T0=300,T,K0,KE;
  double k_B=1.98723e-3;
  double UNITT=418.4070;
  double NfKT,KT;
  double zeta=0.0,V_zeta=0.0,Q,tau=0.01,tau2;
  double PEv,KEv;

  double F,HMFD;

  double pi;

  double *X,*dX_dt,*dF_dX,*mass,massZ;

  double alpha,gamma;

  char *mdinfofilename,*frcfilename,*trjfilename,*velfilename;
  FILE *mdinfofile,*frcfile,*trjfile,*velfile;

  char *line;
  size_t len=0;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  int opt_idx=1;
  
  pi=acos(-1.0);
  
  struct option long_opt[] = {
    {"alpha",1,NULL,'l'},
    {"gamma",1,NULL,'g'},
    {"mass",1,NULL,'m'},
    {"T",1,NULL,'t'},
    {"H",1,NULL,'H'},
    {"tau",1,NULL,'a'},
    {"dt",1,NULL,'x'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };
  
  while((c=getopt_long(argc,argv,"hl:g:m:t:H:a:x:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'l':
      alpha=atof(optarg);
      break;
    case 'g':
      gamma=atof(optarg);
      break;
    case 'm':
      massZ=atof(optarg);
      break;
    case 't':
      T0=atof(optarg);
      break;
    case 'H':
      HMFD=atof(optarg);
      break;
    case 'a':
      tau=atof(optarg);
      break;
    case 'x':
      dt=atof(optarg);
      break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);  exit(1);
    }
  }
  
  progname=*argv;  argc-=optind;  argv+=optind;
  
  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  N = atoi(*argv);
  mdinfofilename= *++argv;
  frcfilename= *++argv;
  trjfilename= *++argv;
  velfilename= *++argv;

  X=(double *)gcemalloc(sizeof(double)*N);
  dX_dt=(double *)gcemalloc(sizeof(double)*N);
  dF_dX=(double *)gcemalloc(sizeof(double)*N);
  mass=(double *)gcemalloc(sizeof(double)*N);

  for (i=0;i<N;++i) mass[i]=massZ;

  mdinfofile=efopen(mdinfofilename,"r");
  frcfile=efopen(frcfilename,"r");

  for (i=0;i<N;++i) {
    fscanf(mdinfofile,"%lf %lf",&X[i],&dX_dt[i]);
    fscanf(frcfile,"%lf",&dF_dX[i]);
  }
  fscanf(mdinfofile,"%lf %lf",&zeta,&V_zeta);

  fclose(mdinfofile);
  fclose(frcfile);

  K0=0.0; for (i=0;i<N;++i) K0+=0.5*mass[i]*dX_dt[i]*dX_dt[i];
  T=K0/(N*k_B)*2.0/UNITT;

  tau=tau/2.0/pi;         
  tau2=tau*tau;          
  KT=k_B*T0;
  NfKT=(N+1)*KT*UNITT;
  Q=tau2*KT*UNITT*N;

  MD_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);

  KE=0.0; for (i=0;i<N;++i) KE+=0.5*mass[i]*dX_dt[i]*dX_dt[i];

  KEv=0.5*Q*V_zeta*V_zeta;
  PEv=NfKT*zeta;

  F=LogMFD_cF(alpha, gamma, HMFD, NfKT, KE,KEv,PEv);

  LogMFD_Propagetor_NH_MP1998_2(X,dX_dt,mass,
				&zeta,&V_zeta,Q,
				NfKT,N,&KE,&KEv,&PEv,
				dt,dt2,nc,wdt4,wdt2,
				&dF_dX,F,alpha,gamma);

  mdinfofile=efopen(mdinfofilename,"w");
  trjfile=efopen(trjfilename,"a");
  velfile=efopen(velfilename,"a");

  for (i=0;i<N;++i) {
    fprintf(trjfile,"%12.8e \n",X[i]);
    fprintf(velfile,"%12.8e \n",dX_dt[i]);
    fprintf(mdinfofile,"%12.8e  %12.8e \n",X[i],dX_dt[i]);
  }
  fprintf(mdinfofile,"%12.8e  %12.8e \n",zeta,V_zeta);

  fclose(mdinfofile);
  fclose(trjfile);
  fclose(velfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename \n",progname);
}


