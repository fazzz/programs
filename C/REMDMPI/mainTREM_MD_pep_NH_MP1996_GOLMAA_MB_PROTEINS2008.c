
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

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"
#include "GOLMAA_MB_PROTEINS2008.h"

#include "REMDMPI.h"
#include "TREM_MD_pep_NH_MP1996_GOLMAA_MB_PROTEINS2008.h"

#include "netcdf_mineL.h"

#define BLANK 0
#define INPF  1
#define TEMP  2
#define TRJ   3
#define OUT   4

#define ON  0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,a;

  int numEX=1,numRE=2;   // REMD

  int vMode=OFF,Equflag=OFF,numstepequ=0;  // flags

  int numatom,numheavyatom,numres,numstep=10000,interval=100;
  double dt=0.001,dt2,wdt2[3],wdt4[3];

  double ep=ep_natatt_hybrid;
  double de=1.0,d=1.0,d2;

  int NCmode=3,nibnum=3,criteria=6.5;

  double pi;

  int nc=1;                                // NHC
  double *T0,*T;                           // NHC
  double k_B=1.98723e-3;                   // NHC
  double UNITT=418.4070;                   // NHC
  double *NfKT,*KT;                        // NHC
  double *zeta,*V_zeta,*Q,tau=0.01,tau2;   // NHC

  double **crd,*refcrd1,*refcrd2,*refcrdAA,*mass,**vel;

  double summass,COM[3];

  double *avePE,*varPE,*aveKE,*varKE,*aveT,*varT;

  struct potential *e;
  struct force *f;
  struct potential_GOLMAA_MB_PROTEINS2008 *e_GOLM;
  //  double p_t=0.0,E_t;

   struct my_netcdf_out_id_MCD *nc_id_MCD,*nc_id_MCD_Equ;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *crdfilename[1000];

  char *outputfilenamebase;
  char *trjfilenamebase,trjfilename[1000],outputfilename[1000],*logfilename,logf[1000];
  char trjfilenameEqu[1000],outputfilenameEqu[1000];

  char *inputfilename,*refcrdfilename1,*refcrdfilename2,*velfilename,*parmfilename;
  char *outputfilename2,*rstfilename="rstcrd",*rstvelfilename="rstvel";

  FILE *inputfile,*crdfile,*velfile,*parmfile,*refcrdfile1,*refcrdfile2;
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
    {"ep",1,NULL,'p'},
    {"cutoff",1,NULL,'c'},
    {"de",1,NULL,'d'},
    {"d",1,NULL,'2'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hv:E:s:e:r::a:i:x:p:c:d:2:",long_opt,&opt_idx))!=-1) {
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
    case 'p':
      ep=atof(optarg);
      break;
    case 'c':
      criteria=atof(optarg);
      break;
    case 'd':
      d=atof(optarg);
      break;
    case '2':
      de=atof(optarg);
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
  //  sprintf(logfilename,"%s.log",progname);

  argc-=optind;
  argv+=optind;

  if (argc < 7) {
    USAGE(progname);
    exit(1);
  }
  inputfilename      = *argv;
  refcrdfilename1    = *++argv;
  refcrdfilename2    = *++argv;
  parmfilename       = *++argv;
  outputfilenamebase = *++argv;
  trjfilenamebase    = *++argv;
  logfilename        = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  j=0;
  for (i=0;i<numatom;++i) {
    if (strncmp(AP.IGRAPH[i],"H",1)==0) {
      ++j;
    }
  }
  numheavyatom=numatom-j;
  numres=AP.NRES;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) {
    //    if (strncmp(AP.IGRAPH[i],"H",1)!=0)
    mass[i]=AP.AMASS[i];
    //    else
    //      mass[i]=0.0;
  }

  refcrd1=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd2=(double *)gcemalloc(sizeof(double)*numatom*3);

  refcrdfile1=efopen(refcrdfilename1,"r");
  getline(&line,&len,refcrdfile1);
  fscanf(refcrdfile1,"%d",&a);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(refcrdfile1,"%lf",&refcrd1[i*3+j]);
  fclose(refcrdfile1);

  refcrdfile2=efopen(refcrdfilename2,"r");
  getline(&line,&len,refcrdfile2);
  fscanf(refcrdfile2,"%d",&a);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(refcrdfile2,"%lf",&refcrd2[i*3+j]);
  fclose(refcrdfile2);

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
  readInputs_TREM_GOLMAA(inputfile,crd,vel,numatom,T0,numRE);
  fclose(inputfile);

  /***************************************************/
  /* for (i=0;i<numRE;++i) {			     */
  /*   for (j=0;j<numatom;++j) {		     */
  /*     if (strncmp(AP.IGRAPH[j],"H",1)==0) 	     */
  /* 	for (k=0;k<3;++k)			     */
  /* 	  crd[i][j*3+k]=0.0;			     */
  /*   }					     */
  /* }						     */
  /***************************************************/

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

  e_GOLM=(struct potential_GOLMAA_MB_PROTEINS2008 *)gcemalloc(sizeof(struct potential_GOLMAA_MB_PROTEINS2008)*numRE);
  
  if (numRE!=num_procs /*|| numRE%2==1*/ ) {
    printf("condition error\n");
    exit(1);
  }

  sprintf(outputfilename,"%s_%d",outputfilenamebase,my_rank+1);
  sprintf(trjfilename,"%s_%d",trjfilenamebase,my_rank+1);

  sprintf(logf,"%s_%d_ex.log",logfilename,my_rank+1);

  tau=tau/2.0/pi;
  tau2=tau*tau;

  KT[my_rank]=k_B*T0[my_rank];
  NfKT[my_rank]=(3.0*numheavyatom+1)*KT[my_rank]*UNITT;
  Q[my_rank]=(tau2)*(KT[my_rank])*UNITT*(3.0*numheavyatom);
  
  ffL_set_calcffandforce(&e[my_rank],&f[my_rank]);
  GOLMAA_MB_PROTEINS2008_ff_calcff_set(&e_GOLM[my_rank],refcrd1,refcrd2,numatom,numres,e[my_rank].parm.indexnb,e[my_rank].parm.numnb,ep,nibnum,criteria);

  d2=d*d;
  d2=d2*k_B;
  de=de*k_B;
  
  MD_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);

  myncL_create_def_MCD(trjfilename,numatom,&(nc_id_MCD[my_rank]));
  outputfile[my_rank]=efopen(outputfilename,"w");

  logfile=efopen(logf,"w");

  if (Equflag==ON) {
    sprintf(outputfilenameEqu,"%s_equ_%d",outputfilenamebase,my_rank+1);
    sprintf(trjfilenameEqu,"%s_equ_%d",trjfilenamebase,my_rank+1);

    myncL_create_def_MCD(trjfilenameEqu,numatom,&(nc_id_MCD_Equ[my_rank]));
    outputfile_Equ[my_rank]=efopen(outputfilenameEqu,"w");

    runMD_pep_NH_MP1996_GOLMAA_MB_PROTEINS2008(crd[my_rank],vel[my_rank],mass,numatom,numheavyatom,
					       &zeta[my_rank],&V_zeta[my_rank],Q[my_rank],
					       e_GOLM[my_rank],e[my_rank],de,d2,
					       T[my_rank],NfKT[my_rank],numstepequ,interval,&l,
					       dt,dt2,wdt2,wdt4,nc,
					       &(avePE[my_rank]),&(aveKE[my_rank]),&(aveT[my_rank]),
					       &(varPE[my_rank]),&(varKE[my_rank]),&(varT[my_rank]),
					       UNITT,k_B,
					       nc_id_MCD_Equ[my_rank],outputfile_Equ[my_rank]);

    fclose(outputfile_Equ[my_rank]);
    nc_close(((nc_id_MCD_Equ[my_rank]).ncid));
  }

  MPI_TREMD_pep_NH_MP1996_GOLMAA_MB_PROTEINS2008(my_rank, num_procs,tag,&status,
						 numEX,numRE,
						 crd,vel,mass,numatom,numheavyatom,
						 zeta,V_zeta,Q,
						 e_GOLM,e,de,d2,
						 T0,NfKT,numstep,interval,
						 dt,dt2,wdt2,wdt4,nc,
						 avePE,aveKE,aveT,
						 varPE,varKE,varT,UNITT,k_B,tau,pi,
						 nc_id_MCD,outputfile,logfile);

  fclose(logfile);
  
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

