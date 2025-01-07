
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "mpi.h"

#include "PTL.h"
#include "EF.h"
#include "MDrun.h"
#include "MDrun.h"
#include "MD_NHC_MP1996.h"

#include "REMDMPI.h"

#include "netcdf_mineL.h"

#define ON  0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0;

  int numEX=1,numRE=2;   // REMD

  int vMode=OFF,Equflag=OFF,numstepequ=0;  // flags

  int numatom,numstep=10000,interval=100;
  double dt=0.001,dt2,wdt2[3],wdt4[3];

  double pi;

  int nc=1;                                // NHC
  double *T0,*T;                           // NHC
  double k_B=1.98723e-3;                   // NHC
  double UNITT=418.4070;                   // NHC
  double *NfKT,*KT;                        // NHC
  double *zeta,*V_zeta,*Q,tau=0.01,tau2;   // NHC

  double **crd,**vel,*mass;

  double *avePE,*varPE,*aveKE,*varKE,*aveT,*varT;

  struct potential *e;
  struct force *f;

  struct my_netcdf_out_id_MCD *nc_id_MCD,*nc_id_MCD_Equ;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename;
  char *crdfilename[1000],*velfilename,*parmfilename;

  char *outputfilenamebase;
  char *trjfilenamebase,trjfilename[1000],outputfilename[1000],logfilename[100];
  char trjfilenameEqu[1000],outputfilenameEqu[1000];

  FILE *inputfile,*crdfile,*velfile,*parmfile;
  FILE **outputfile,**outputfile_Equ,*logfile;

  char *progname;

  int opt_idx=1;

  int my_rank,num_procs,tag = 0;  // MPI
  MPI_Status status;              // MPI

  // MPI
  MPI_Init(&argc, &argv);
  
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  // MPI

  pi=acos(-1.0);

  struct option long_opt[] = {
    {"vMode",1,NULL,'v'},
    {"equ",1,NULL,'E'},
    {"nums",1,NULL,'s'},
    {"numEX",1,NULL,'e'},
    {"numRE",1,NULL,'r'},
    {"tau",1,NULL,'a'},
    {"int",1,NULL,'i'},
    {"dt",1,NULL,'x'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hv:E:s:e:r::a:i:x:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 's':
      numstep=atoi(optarg);
      break;
    case 'E':
      Equflag=ON;
      numstepequ=atoi(optarg);
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
    case 'e':
      numEX=atoi(optarg);
      break;
    case 'r':
      numRE=atoi(optarg);
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  inputfilename      = *argv;
  parmfilename       = *++argv;
  outputfilenamebase = *++argv;
  trjfilenamebase    = *++argv;

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
  T=(double *)gcemalloc(sizeof(double)*numRE);
  KT=(double *)gcemalloc(sizeof(double)*numRE);
  NfKT=(double *)gcemalloc(sizeof(double)*numRE);
  Q=(double *)gcemalloc(sizeof(double)*numRE);
  nc_id_MCD=(struct my_netcdf_out_id_MCD *)gcemalloc(sizeof(struct my_netcdf_out_id_MCD)*numRE);
  nc_id_MCD_Equ=(struct my_netcdf_out_id_MCD *)gcemalloc(sizeof(struct my_netcdf_out_id_MCD)*numRE);
  outputfile=(FILE **)gcemalloc(sizeof(FILE *)*numRE);
  outputfile_Equ=(FILE **)gcemalloc(sizeof(FILE *)*numRE);

  avePE=(double *)gcemalloc(sizeof(double)*numRE);
  aveKE=(double *)gcemalloc(sizeof(double)*numRE);
  aveT=(double *)gcemalloc(sizeof(double)*numRE);

  varPE=(double *)gcemalloc(sizeof(double)*numRE);
  varKE=(double *)gcemalloc(sizeof(double)*numRE);
  varT=(double *)gcemalloc(sizeof(double)*numRE);

  inputfile=efopen(inputfilename,"r");
  readInputs(inputfile,crd,vel,numatom,T0);
  fclose(inputfile);

  outputfile=(FILE **)gcemalloc(sizeof(FILE *));

  zeta=(double *)gcemalloc(sizeof(double)*numRE);
  V_zeta=(double *)gcemalloc(sizeof(double)*numRE);

  if ( vMode==OFF ) {
    for (i=0;i<numRE;++i) {
      MD_Generate_inivelo(vel[i],mass,numatom,k_B*T0[i]*UNITT);

      zeta[i]=0.0;
      V_zeta[i]=0.0;
    }
  }
  
  e=(struct potential *)gcemalloc(sizeof(struct potential)*numRE);
  f=(struct force *)gcemalloc(sizeof(struct force)*numRE);

  if (numRE!=num_procs || numRE%2==1 ) {
    printf("condition error\n");
    exit(1);
  }
  
  sprintf(outputfilename,"%s_%d",outputfilenamebase,my_rank+1);
  sprintf(trjfilename,"%s_%d",trjfilenamebase,my_rank+1);

  tau=tau/2.0/pi;
  tau2=tau*tau;

  KT[my_rank]=k_B*T0[my_rank];
  NfKT[my_rank]=(3.0*numatom+1)*KT[my_rank]*UNITT;
  Q[my_rank]=(tau2)*(KT[my_rank])*UNITT*(3.0*numatom);
  
  ffL_set_calcffandforce(&e[my_rank],&f[my_rank]);
  
  MD_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);

  myncL_create_def_MCD(trjfilename,numatom,&(nc_id_MCD[my_rank]));
  outputfile[my_rank]=efopen(outputfilename,"w");


  if (Equflag==ON) {
    sprintf(outputfilenameEqu,"%s_equ_%d",outputfilenamebase,my_rank+1);
    sprintf(trjfilenameEqu,"%s_qu_%d",trjfilenamebase,my_rank+1);

    myncL_create_def_MCD(trjfilenameEqu,numatom,&(nc_id_MCD_Equ[my_rank]));
    outputfile_Equ[my_rank]=efopen(outputfilenameEqu,"w");

    runMD_NHC_MP1998_Amber_AAFF(crd[my_rank],vel[my_rank],mass,numatom,
				&zeta[my_rank],&V_zeta[my_rank],Q[my_rank],
				e[my_rank],f[my_rank],
				T[my_rank],NfKT[my_rank],numstepequ,interval,&l,
				dt,dt2,wdt2,wdt4,nc,
				&(avePE[my_rank]),&(aveKE[my_rank]),&(aveT[my_rank]),
				&(varPE[my_rank]),&(varKE[my_rank]),&(varT[my_rank]),
				UNITT,k_B,
				nc_id_MCD_Equ[my_rank],outputfile_Equ[my_rank]);

    fclose(outputfile_Equ[my_rank]);
    nc_close(((nc_id_MCD_Equ[my_rank]).ncid));
  }

  /***************************************************************************************************/
  /* MD_Init(numatom,T0[my_rank],k_B,UNITT,tau,&tau2,&KT[my_rank],&NfKT[my_rank],&Q[my_rank],	     */
  /* 	  &e[my_rank],&f[my_rank],nc,dt,&dt2,wdt2,wdt4,pi,					     */
  /* 	  nc_id_MCD[my_rank],trjfilename,outputfile[my_rank],outputfilename);			     */
  /***************************************************************************************************/
  
  MPI_TREMD_pep_NHC_MP1998_Amber_AAFF(my_rank, num_procs,tag,&status,
				      numEX,numRE,
				      crd,vel,mass,numatom,
				      zeta,V_zeta,Q,e,f,
				      T0,NfKT,numstep,interval,
				      dt,dt2,wdt2,wdt4,nc,
				      avePE,aveKE,aveT,
				      varPE,varKE,varT,UNITT,k_B,tau,pi,
				      nc_id_MCD,outputfile);
  
  /************************************************************************/
  /* MD_Fine(nc_id_MCD[my_rank],outputfile[my_rank],logfilename,	  */
  /* 	  &(avePE[my_rank]),&(aveKE[my_rank]),&(aveT[my_rank]),		  */
  /* 	  &(varPE[my_rank]),&(varKE[my_rank]),&(varT[my_rank]),		  */
  /* 	  UNITT,k_B);							  */
  /************************************************************************/

  fclose(outputfile[my_rank]);
  nc_close(((nc_id_MCD[my_rank]).ncid));
  
  MPI_Finalize();
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parmfilename outputfilenamebase trjfilenamebase\n",progname);
}
