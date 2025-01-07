
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

  ///////////////// TACCM //////////////////////
  int massflag=OFF;
  double massX=1.0;
  ///////////////// TACCM //////////////////////

  int nc=1;                                // NHC
  double T0,*T;                           // NHC
  double k_B=1.98723e-3;                   // NHC
  double UNITT=418.4070;                   // NHC
  double NfKT,KT;                        // NHC
  double *zeta,*V_zeta,Q,tau=0.01,tau2;   // NHC

  double **crd,**vel,*mass;

  double *avePE,*varPE,*aveKE,*varKE,*aveT,*varT;

  struct potential *e;
  struct force *f;

  struct my_netcdf_out_id_MCD *nc_id_MCD,*nc_id_MCD_Equ;

  ///////////////// TACCM //////////////////////
  double **theta;
  double **Z,**velZ,**accZ;
  int numZ;
  double *TobjZ,*KEobjZ,*KBTZ,*TZ;
  double massZ=100.0;
  double **frcZ;
  double *KEZ,*PEZ,*KEvZ,*PEvZ,*EtZ;
  double *zetaZ,*V_zetaZ,*QZ,*NfKTZ;
  double KZ=10.0;
  int *indexTACCM,**pairsZ;
  char *trjfileZbase,*trjfileThetaZ;
  char *TACCMfilename,trjfilenameZ[1000],trjfilenameTheta[1000];
  FILE *TACCMfile,**trjfileZ,**trjfileTheta;

  double tauZ,tau2Z;

  double *avePEZ,*varPEZ,*aveKEZ,*varKEZ,*aveTZ,*varTZ;
  ///////////////// TACCM //////////////////////

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
    {"temp",1,NULL,'T'},
    {"numEX",1,NULL,'e'},
    {"numRE",1,NULL,'r'},
    {"tau",1,NULL,'a'},
    {"int",1,NULL,'i'},
    {"dt",1,NULL,'x'},
    {"T",1,NULL,'T'},
    // TACCM ///////////////
    {"mZ",1,NULL,'m'},
    {"KZ",1,NULL,'K'},
    {"massX",1,NULL,'X'},
    // TACCM ///////////////
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hv:E:s:T:e:r::a:i:x:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 's':
      numstep=atoi(optarg);
      break;
    case 'E':
      Equflag=ON;
      numstepequ=atoi(optarg);
      break;
    case 'T':
      T0=atof(optarg);
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
      ///////////////// TACCM //////////////////////
    case 'm':
      massZ=atof(optarg);
      break;
    case 'K':
      KZ=atof(optarg);
      break;
    /*************************/
    /* case 'B':	     */
    /*   TobjZ=atof(optarg); */
    /*   break;		     */
    /*************************/
    case 'X':
      massflag=ON;
      massX=atof(optarg);
      break;
      ///////////////// TACCM //////////////////////
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

  if (argc < 7) {
    USAGE(progname);
    exit(1);
  }
  inputfilename      = *argv;
  parmfilename       = *++argv;
  ///////////////// TACCM //////////////////////
  TACCMfilename  = *++argv;
  ///////////////// TACCM //////////////////////
  outputfilenamebase = *++argv;
  trjfilenamebase    = *++argv;
  ///////////////// TACCM //////////////////////
  trjfileZbase       = *++argv;
  trjfileThetaZ      = *++argv;
  ///////////////// TACCM //////////////////////

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
  //  T0=(double *)gcemalloc(sizeof(double)*numRE);
  T=(double *)gcemalloc(sizeof(double)*numRE);
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

  ///////////////// TACCM //////////////////////////
  TobjZ=(double *)gcemalloc(sizeof(double)*numRE);
  ///////////////// TACCM //////////////////////////

  inputfile=efopen(inputfilename,"r");
  readInputs(inputfile,crd,vel,numatom,TobjZ);
  fclose(inputfile);

  outputfile=(FILE **)gcemalloc(sizeof(FILE *));

  zeta=(double *)gcemalloc(sizeof(double)*numRE);
  V_zeta=(double *)gcemalloc(sizeof(double)*numRE);

  if ( vMode==OFF ) {
    for (i=0;i<numRE;++i) {
      MD_Generate_inivelo(vel[i],mass,numatom,k_B*T0*UNITT);

      zeta[i]=0.0;
      V_zeta[i]=0.0;
    }
  }
  
  e=(struct potential *)gcemalloc(sizeof(struct potential)*numRE);
  f=(struct force *)gcemalloc(sizeof(struct force)*numRE);


  ///////////////// TACCM //////////////////////
  TACCMfile=efopen(TACCMfilename,"r");
  fscanf(TACCMfile,"%d",&numZ);
  pairsZ=(int **)gcemalloc(sizeof(int *)*numZ);
  for (i=0;i<numZ;++i) pairsZ[i]=(int *)gcemalloc(sizeof(int)*5);
  for (i=0;i<numZ;++i) {
    for (j=0;j<4;++j) 
      fscanf(TACCMfile,"%d",&pairsZ[i][j]);
    fscanf(TACCMfile,"%d",&pairsZ[i][j]);
  }
  fclose(TACCMfile);
  theta=(double **)gcemalloc(sizeof(double *)*numRE);
  Z=(double **)gcemalloc(sizeof(double *)*numRE);
  velZ=(double **)gcemalloc(sizeof(double *)*numRE);
  frcZ=(double **)gcemalloc(sizeof(double)*numRE);
  for (i=0;i<numRE;++i) {
    theta[i]=(double *)gcemalloc(sizeof(double)*numZ);
    Z[i]=(double *)gcemalloc(sizeof(double)*numZ);
    velZ[i]=(double *)gcemalloc(sizeof(double)*numZ);
    frcZ[i]=(double *)gcemalloc(sizeof(double)*numZ);
  }
  TACCM_CTheta(crd[my_rank],numatom,theta[my_rank],numZ,pairsZ,pi);
  for (j=0;j<numZ;++j) Z[my_rank][j]=theta[my_rank][j];
  /************************************************/
  /* for (j=0;j<numZ;++j) 			  */
  /*   printf("%d %8.4lf\n",j,theta[my_rank][j]); */
  /************************************************/
  KBTZ=(double *)gcemalloc(sizeof(double)*numRE); 
  NfKTZ=(double *)gcemalloc(sizeof(double)*numRE);
  QZ=(double *)gcemalloc(sizeof(double)*numRE);
  KEZ=(double *)gcemalloc(sizeof(double)*numRE);
  PEZ=(double *)gcemalloc(sizeof(double)*numRE);
  EtZ=(double *)gcemalloc(sizeof(double)*numRE);
  zetaZ=(double *)gcemalloc(sizeof(double)*numRE);
  V_zetaZ=(double *)gcemalloc(sizeof(double)*numRE);
  trjfileZ=(FILE **)gcemalloc(sizeof(FILE *)*numRE);
  trjfileTheta=(FILE **)gcemalloc(sizeof(FILE *)*numRE);
  for (i=0;i<numRE;++i) {
    zetaZ[i]=0.0;
    V_zetaZ[i]=0.0;
    KBTZ[i]=k_B*TobjZ[i];
    NfKTZ[i]=(numZ+1)*KBTZ[i]*UNITT;
    KEZ[i]=TACCM_MD_Generate_inivelo(velZ[i],massZ,numZ,k_B*TobjZ[i]*UNITT);
  }
  avePEZ=(double *)gcemalloc(sizeof(double)*numRE);
  varPEZ=(double *)gcemalloc(sizeof(double)*numRE);
  aveKEZ=(double *)gcemalloc(sizeof(double)*numRE);
  varKEZ=(double *)gcemalloc(sizeof(double)*numRE);
  aveTZ=(double *)gcemalloc(sizeof(double)*numRE);
  varTZ=(double *)gcemalloc(sizeof(double)*numRE);
  ///////////////// TACCM //////////////////////


  if (numRE!=num_procs || numRE%2==1 ) {
    printf("condition error\n");
    exit(1);
  }
  
  sprintf(outputfilename,"%s_%d",outputfilenamebase,my_rank+1);
  sprintf(trjfilename,"%s_%d",trjfilenamebase,my_rank+1);
  ///////////////// TACCM //////////////////////
  sprintf(trjfilenameZ,"%s_%d",trjfileZbase,my_rank+1);
  sprintf(trjfilenameTheta,"%s_%d",trjfileThetaZ,my_rank+1);
  ///////////////// TACCM //////////////////////

  tau=tau/2.0/pi;
  tau2=tau*tau;
  
  ///////////////// TACCM //////////////////////
  for (i=0;i<numRE;++i) {
    QZ[i]=tau2*KBTZ[i]*UNITT*numZ;
  }
  ///////////////// TACCM //////////////////////

  KT=k_B*T0;
  NfKT=(3.0*numatom+1)*KT*UNITT;
  Q=tau2*KT*UNITT*(3.0*numatom);
  
  ffL_set_calcffandforce(&e[my_rank],&f[my_rank]);
  
  MD_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);

  myncL_create_def_MCD(trjfilename,numatom,&(nc_id_MCD[my_rank]));
  outputfile[my_rank]=efopen(outputfilename,"w");
  ///////////////// TACCM //////////////////////
  trjfileZ[my_rank]=efopen(trjfilenameZ,"w");
  trjfileTheta[my_rank]=efopen(trjfilenameTheta,"w");
  ///////////////// TACCM //////////////////////

  if (Equflag==ON) {
    ;
  }

  /***************************************************************************************************/
  /* MD_Init(numatom,T0[my_rank],k_B,UNITT,tau,&tau2,&KT[my_rank],&NfKT[my_rank],&Q[my_rank],	     */
  /* 	  &e[my_rank],&f[my_rank],nc,dt,&dt2,wdt2,wdt4,pi,					     */
  /* 	  nc_id_MCD[my_rank],trjfilename,outputfile[my_rank],outputfilename);			     */
  /***************************************************************************************************/

  /****************************************************************************************************/
  /* printf("yes This is line 360 in main%d\n",my_rank);					      */
  /* 												      */
  /* printf("Kapa=%8.4lf massZ=%8.4lf QZ=%8.4lf NfKTZ=%8.4lf\n",KZ,massZ,QZ[my_rank],NfKTZ[my_rank]); */
  /****************************************************************************************************/
  
  MPI_TREM_TACCM_MD_pep_NHC_MP1998_Amber_AAFF(my_rank, num_procs,tag,&status,
  					      numEX,numRE,
  					      crd,vel,mass,numatom,
  					      zeta,V_zeta,Q,e,f,
  					      T0,NfKT,numstep,interval,
  					      dt,dt2,wdt2,wdt4,nc,
  					      avePE,aveKE,aveT,
  					      varPE,varKE,varT,
  					      UNITT,k_B,tau,pi,
  					      nc_id_MCD,outputfile,
					      //////////////// TACCM ///////////////////////
  					      Z,velZ,massZ,
  					      zetaZ,V_zetaZ,
  					      TobjZ,QZ,NfKTZ,numZ,
  					      KZ,pairsZ,
  					      avePEZ,aveKEZ,aveTZ,
  					      varPEZ,varKEZ,varTZ,
  					      trjfileZ,trjfileTheta
					      //////////////// TACCM ///////////////////////
					      );
  
  /************************************************************************/
  /* MD_Fine(nc_id_MCD[my_rank],outputfile[my_rank],logfilename,	  */
  /* 	  &(avePE[my_rank]),&(aveKE[my_rank]),&(aveT[my_rank]),		  */
  /* 	  &(varPE[my_rank]),&(varKE[my_rank]),&(varT[my_rank]),		  */
  /* 	  UNITT,k_B);							  */
  /************************************************************************/

  fclose(outputfile[my_rank]);
  nc_close(((nc_id_MCD[my_rank]).ncid));
  fclose(trjfileZ[my_rank]);
  fclose(trjfileTheta[my_rank]);
  
  MPI_Finalize();
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parmfilename outputfilenamebase trjfilenamebase\n",progname);
}
