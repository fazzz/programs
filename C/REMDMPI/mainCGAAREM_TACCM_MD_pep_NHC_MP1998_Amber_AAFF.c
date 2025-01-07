
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "mpi.h"

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

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
  int i,j,k,l=0,d;
  int numEX=1,numRE=2;
  int vMode=OFF,Equflag=OFF,numstepequ=0;

  int numatom,numheavyatom,numres,numstep=10000,interval=100;
  double dt=0.001,dt2,wdt2[3],wdt4[3];

  double pi;

  ///////////////// TACCM //////////////////////
  int massflag=OFF;
  double massX=1.0;
  ///////////////// TACCM //////////////////////

  int nc=1;                          
  double T0AA,T0CG;
  double k_B=1.98723e-3;             
  double UNITT=418.4070;             
  double KTAA,KTCG,tau=0.01,tau2;                    

  double TAA,TCG;                       

  double zetaAA,V_zetaAA,QAA,NfKTAA;
  double zetaCG,V_zetaCG,QCG,NfKTCG;

  double *crdAA,*velAA,*mass;
  double *crdCG,*velCG,*refcrdCG;

  double avePEAA,varPEAA,aveKEAA,varKEAA,aveTAA,varTAA;
  double avePECG,varPECG,aveKECG,varKECG,aveTCG,varTCG;

  double ep=ep_natatt_hybrid;

  int NCmode=3,nibnum=3,criteria=6.5;

  struct potential e;
  struct force f;
  struct potential_GOLMAA_PROTEINS2008 e_GOLM;
  int numnb,num14;

  struct my_netcdf_out_id_MCD nc_id_MCDAA;
  struct my_netcdf_out_id_MCD nc_id_MCDCG;

  ///////////////// TACCM //////////////////////
  double **theta;
  double **Z,**velZ,**accZ;
  int numZ;
  double TobjZ,*KEobjZ,KBTZ,*TZ;
  double massZ=100.0;
  double **frcZ;
  double *KEZ,*PEZ,*KEvZ,*PEvZ,*EtZ;
  double *zetaZ,*V_zetaZ,QZ,NfKTZ;
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
  char *refcrdfilename;

  char *outputfilenameAA;
  char *trjfilenameAA;
  char *outputfilenameCG;
  char *trjfilenameCG;

  FILE *inputfile,*crdfile,*refcrdfile,*velfile,*parmfile,*logfilename;

  FILE *outputfileAA,*logfileAA;
  FILE *outputfileCG,*logfileCG;

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
    {"tau",1,NULL,'a'},
    {"int",1,NULL,'i'},
    {"dt",1,NULL,'x'},
    {"TAA",1,NULL,'T'},
    {"TCG",1,NULL,'t'},
    {"TB",1,NULL,'B'},
    // TACCM ///////////////
    {"mZ",1,NULL,'m'},
    {"KZ",1,NULL,'K'},
    {"massX",1,NULL,'X'},
    // TACCM ///////////////
    // GOLMAA //////////////
    {"ep",1,NULL,'p'},
    {"cutoff",1,NULL,'c'},
    // GOLMAA //////////////
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hv:E:s:e:a:i:x:T:t:B:m:K:X:p:c:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 's':
      numstep=atoi(optarg);
      break;
    case 'E':
      Equflag=ON;
      numstepequ=atoi(optarg);
      break;
    case 'T':
      T0AA=atof(optarg);
      break;
    case 't':
      T0CG=atof(optarg);
      break;
    case 'B':
      TobjZ=atof(optarg);
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
      ///////////////// TACCM //////////////////////
    case 'm':
      massZ=atof(optarg);
      break;
    case 'K':
      KZ=atof(optarg);
      break;
    case 'X':
      massflag=ON;
      massX=atof(optarg);
      break;
      ///////////////// TACCM //////////////////////
      ///////////////// GOLMAA //////////////////////
    case 'p':
      ep=atof(optarg);
      break;
    case 'c':
      criteria=atof(optarg);
      break;
      ///////////////// GOLMAA //////////////////////
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

  if (argc < 10) {
    USAGE(progname);
    exit(1);
  }
  inputfilename      = *argv;
  ///////////////// GOLMAA /////////////////////
  refcrdfilename    = *++argv;
  ///////////////// GOLMAA /////////////////////
  parmfilename       = *++argv;
  ///////////////// TACCM //////////////////////
  TACCMfilename  = *++argv;
  ///////////////// TACCM //////////////////////
  outputfilenameAA = *++argv;
  trjfilenameAA    = *++argv;
  outputfilenameCG = *++argv;
  trjfilenameCG    = *++argv;
  ///////////////// TACCM //////////////////////
  trjfileZbase       = *++argv;
  trjfileThetaZ      = *++argv;
  ///////////////// TACCM //////////////////////

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  ///////////////// GOLMAA /////////////////////
  j=0;
  for (i=0;i<numatom;++i) {
    if (strncmp(AP.IGRAPH[i],"H",1)==0) {
      ++j;
    }
  }
  numheavyatom=numatom-j;
  numres=AP.NRES;
  ///////////////// GOLMAA /////////////////////
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];
  
  crdAA=(double *)gcemalloc(sizeof(double)*numatom*3);
  velAA=(double *)gcemalloc(sizeof(double)*numatom*3);
  ///////////////// GOLMAA /////////////////////
  crdCG=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrdCG=(double *)gcemalloc(sizeof(double)*numatom*3);
  velCG=(double *)gcemalloc(sizeof(double)*numatom*3);
  ///////////////// GOLMAA /////////////////////

  ///////////////// TACCM //////////////////////////
  //  TobjZ=(double *)gcemalloc(sizeof(double)*numRE);
  ///////////////// TACCM //////////////////////////

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) { 
    for (j=0;j<3;++j) { 
      fscanf(inputfile,"%lf",&crdAA[i*3+j]);
      crdCG[i*3+j]=crdAA[i*3+j];
    }
  }
  fclose(inputfile);

  refcrdfile=efopen(refcrdfilename,"r");
  getline(&line,&len,refcrdfile);
  fscanf(refcrdfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(refcrdfile,"%lf",&refcrdCG[i*3+j]);
  fclose(refcrdfile);

  if ( vMode==OFF ) {
    MD_Generate_inivelo(velAA,mass,numatom,k_B*T0AA*UNITT);
    MD_Generate_inivelo(velCG,mass,numatom,k_B*T0CG*UNITT);
    for (i=0;i<numatom;++i) 
      if (strncmp(AP.IGRAPH[i],"H",1)==0) 
	for (j=0;j<3;++j) 
	  velCG[i*3+j]=0.0;

    zetaAA=0.0;
    V_zetaAA=0.0;
    zetaCG=0.0;
    V_zetaCG=0.0;
  }
  
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
  TACCM_CTheta(crdAA,numatom,theta[0],numZ,pairsZ,pi);
  TACCM_CTheta(crdCG,numatom,theta[1],numZ,pairsZ,pi);
  for (j=0;j<numZ;++j) Z[my_rank][j]=theta[my_rank][j];

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
    KEZ[i]=TACCM_MD_Generate_inivelo(velZ[i],massZ,numZ,k_B*TobjZ*UNITT);
  }
  KBTZ=k_B*TobjZ;
  NfKTZ=(numZ+1)*KBTZ*UNITT;
  avePEZ=(double *)gcemalloc(sizeof(double)*numRE);
  varPEZ=(double *)gcemalloc(sizeof(double)*numRE);
  aveKEZ=(double *)gcemalloc(sizeof(double)*numRE);
  varKEZ=(double *)gcemalloc(sizeof(double)*numRE);
  aveTZ=(double *)gcemalloc(sizeof(double)*numRE);
  varTZ=(double *)gcemalloc(sizeof(double)*numRE);
  ///////////////// TACCM //////////////////////

  if ( num_procs != 2 ) {
    printf("condition error\n");
    exit(1);
  }
  
  ///////////////// TACCM //////////////////////
  sprintf(trjfilenameZ,"%s_%d",trjfileZbase,my_rank+1);
  sprintf(trjfilenameTheta,"%s_%d",trjfileThetaZ,my_rank+1);
  ///////////////// TACCM //////////////////////

  tau=tau/2.0/pi;
  tau2=tau*tau;
  
  ///////////////// TACCM //////////////////////
  QZ=tau2*KBTZ*UNITT*numZ;
  ///////////////// TACCM //////////////////////

  KTAA=k_B*T0AA;
  KTCG=k_B*T0CG;
  NfKTAA=(3.0*numatom+1)*KTAA*UNITT;
  QAA=tau2*KTAA*UNITT*(3.0*numatom);
  NfKTCG=(3.0*numheavyatom+1)*KTCG*UNITT;
  QCG=tau2*KTCG*UNITT*(3.0*numheavyatom);
  
  ffL_set_calcffandforce(&e,&f);

  ffL_set_non_bonding_index_1(&numnb,&num14);
  e.parm.numnb=numnb;
  e.parm.num14=num14;
  e.parm.indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  e.parm.index14=(int *)gcemalloc(sizeof(int)*num14*2);
  ffL_set_non_bonding_index_2(e.parm.indexnb,e.parm.index14);

  GOLMAA_PROTEINS2008_ff_set_calcff_b(&e_GOLM,refcrdCG,numatom,numres,
				      e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);

  ffL_set_calcffandforce(&e,&f);

  MD_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);

  myncL_create_def_MCD(trjfilenameAA,numatom,&nc_id_MCDAA);
  myncL_create_def_MCD(trjfilenameCG,numatom,&nc_id_MCDCG);
  outputfileAA=efopen(outputfilenameAA,"w");
  outputfileCG=efopen(outputfilenameCG,"w");

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

  MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_Amber_AAFF(my_rank, num_procs,tag,&status,numEX,
    						  ////////////// AA ////////////////////////////
    						  crdAA,velAA,mass,numatom, zetaAA,V_zetaAA,QAA,e,f,NfKTAA,
    						  &avePEAA,&aveKEAA,&aveTAA, &varPEAA,&varKEAA,&varTAA,
    						  nc_id_MCDAA,outputfileAA,
    						  ////////////// CG ////////////////////////////
    						  crdCG,velCG, zetaCG,V_zetaCG,QCG,e_GOLM,NfKTCG,
    						  &avePECG,&aveKECG,&aveTCG, &varPECG,&varKECG,&varTCG,
   						  nc_id_MCDCG,outputfileCG,
    						  ////////////// COMMON /////////////////////////
    						  dt,dt2,wdt2,wdt4,nc, numstep,interval,T0AA,T0CG, UNITT,k_B,tau,pi,
    						  //////////////// TACCM ///////////////////////
    						  Z,velZ,massZ, zetaZ,V_zetaZ,
    						  TobjZ,QZ,NfKTZ,numZ,  KZ,pairsZ,
    						  avePEZ,aveKEZ,aveTZ, varPEZ,varKEZ,varTZ,
   						  trjfileZ,trjfileTheta);
  
  /************************************************************************/
  /* MD_Fine(nc_id_MCD[my_rank],outputfile[my_rank],logfilename,	  */
  /* 	  &(avePE[my_rank]),&(aveKE[my_rank]),&(aveT[my_rank]),		  */
  /* 	  &(varPE[my_rank]),&(varKE[my_rank]),&(varT[my_rank]),		  */
  /* 	  UNITT,k_B);							  */
  /************************************************************************/

  fclose(outputfileAA);
  fclose(outputfileCG);
  nc_close(nc_id_MCDAA.ncid);
  nc_close(nc_id_MCDCG.ncid);
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
