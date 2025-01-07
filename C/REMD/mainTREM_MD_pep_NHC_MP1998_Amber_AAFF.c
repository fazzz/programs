
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#include "PT.h"
#include "EF.h"

#include "REMDMPI.h"

#define BLANK 0
#define INPF  1
#define TEMP  2

int main(int argc, char *argv[]) {
  int i,j,k,d;
  double f;

  int numEX=1,numRE=1;

  int State,Tflag=OFF;

  int numatom,numstep=10000,interval=100;
  double dt=0.001,dt2,wdt2[3],wdt4[3];

  double pi;

  int nc=1;
  double T0=300,T,K0,KE;
  double k_B=1.98723e-3;
  double UNITT=418.4070;
  double NfKT,KT;
  double *zeta,*V_zeta,Q,tau=0.01,tau2;
  double *PEv,*KEv;

  double **crd,**vel,*mass;

  double *avePE,*varPE,*aveKE,*varKE,*aveT,*varT;

  struct potential *e;
  struct force *f;

  struct my_netcdf_out_id_MCD *nc_id_MCD;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*crdfilename,*velfilename,*parmfilename;
  char *trjfilename,*outputfilename,logfilename[100];

  FILE *inputfile,*crdfile,*velfile,*parmfile;
  FILE **outputfile,*logfile;

  char *progname;

  int opt_idx=1;

  int my_rank,num_procs,tag = 0;
  MPI_Status status;

  pi=acos(-1.0);

  struct option long_opt[] = {
    {"vMode",1,NULL,'v'},
    {"nums",1,NULL,'s'},
    {"numEX",1,NULL,'e'},
    {"numRE",1,NULL,'r'},
    {"tau",1,NULL,'a'},
    {"int",1,NULL,'i'},
    {"dt",1,NULL,'x'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"*hs:a:i:x:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 's':
      numstep=atoi(optarg);
      break;
    case 'v':
      vMode=ON;
      velfilename=optarg;
      break;
    case 'a':
      tau=atof(optarg);
      break;
    case 'i':
      interval=atoi(optarg);
      break;
    case 'x':
      dt=atof(optarg);
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
  sprintf(logfilename,"%s.log",progname);

  argc-=optind;
  argv+=optind;

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  inputfilename     = *argv;
  parmfilename      = *++argv;
  outputfilename    = *++argv;
  trjfilename       = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 

  numatom=AP.NATOM;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];
  
  crd=(double **)gcemalloc(sizeof(double *)*numRE);
  vel=(double **)gcemalloc(sizeof(double *)*numRE);
  for (i=0;i<numRE;++i) {
    crd[i]=(double *)gcemalloc(sizeof(double)*numatom*3);
    vel[i]=(double *)gcemalloc(sizeof(double)*numatom*3);
  }
  T0=(double *)gcemalloc(sizeof(double)*numRE);
  nc_id_MCD=(struct my_netcdf_out_id_MCD *)gcemalloc(sizeof(struct my_netcdf_out_id_MCD)*numRE);
  outputfile=(FILE **)gcemalloc(sizeof(FILE *)*numRE);

  avePE=(double *)gcemalloc(sizeof(double)*numRE);
  aveKE=(double *)gcemalloc(sizeof(double)*numRE);
  aveT=(double *)gcemalloc(sizeof(double)*numRE);

  varPE=(double *)gcemalloc(sizeof(double)*numRE);
  varKE=(double *)gcemalloc(sizeof(double)*numRE);
  varT=(double *)gcemalloc(sizeof(double)*numRE);


  i=0;  j=0;  f=0.0;
  State=INPF;
  Tflag=OFF;
  inputfilename=efopen(inputfile,"r");
  while ((c=getc(inputfile))!=-1){
    if (State=INPF) {
      if (c==' ' || c=='\n') {
	if (j>0) {
	  crdfilename[j]='\0';
	  crdfile=efopen(crdfilename,"r");
	  getline(&line,&len,crdfile);
	  fscanf(crdfile,"%d",&d);
	  for (k=0;k<numatom*3;++k) fscanf(crdfile,"%lf",&crd[i][k]);
	  fclose(crdfile);
	  State=TEMP;
	  j=0;
	}
      }
      else {
	crdfilename[j]=c;
	++j;
      }
    }
    else if (State==TEMP) {
      if (isdigit(c)) {
	d=(c-'0');
	f=f*10+(double)d;
	Tflag=ON;
      }
      else {
	printf("error\n");
	exit(1);
      }
      if (c==' ' || c=='\n') {
	if (Tflag==ON) {
	  T[i]=f;
	  ++i;
	  State=INPF;
	  Tflag=OFF;
	}
      }
    }
  }    
  fclose(inputfile);

  zeta=(double *)gcemalloc(sizeof(double)*numRE);
  V_zeta=(double *)gcemalloc(sizeof(double)*numRE);

  if ( vMode==OFF ) {
    for (i=0;i<numRE;++i) {
      MD_Generate_inivelo(vel[i],mass,numatom,k_B*T0*UNITT);

      zeta[i]=0.0;
      V_zeta[i]=0.0;
    }
  }

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &p_num);

  if (numRE!=num_procs) {
    printf("condition error\n");
    exit(1);
  }

  MD_Init(numatom,T0[numrank],k_B,UNITT,tau,&tau2,&KT[numrank],&NfKT[numrank],&Q,
	  &e[numrank],&f[numrank],nc,dt,&dt2,wdt2,wdt4,pi);

  TREMD_pep_NHC_MP1998_Amber_AAFF(numrank, num_procs,tag,status,
				  numEX,numRE,
				  crd,vel,mass,numatom,
				  zeta,V_zeta,Q,e,f,
				  T,numstep,interval,
				  dt,dt2,wdt2,wdt4
				  avePE,aveKE,aveT,
				  varPE,varKE,varT,UNITT,k_B,
				  nc_id_MCD,outputfile);

  MD_Fine(nc_id_MCD[numrank],outputfile[numrank],logfilename,
	  avePE[numrank],aveKE[numrank],aveT[numrank],
	  varPE[numrank],varKE[numrank],varT[numrank],
	  UNITT,k_B);

  MPI_Finalize();
  
  return 0;
}
