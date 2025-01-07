
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "REMDCGAA_TACCM_MPI_2_Amber_PROTEINS2008.h"
#include "REMDCGAA_TACCM_MPI_2_Amber_2PROTEINS2008.h"
#include "REMDCGAA_TACCM_MPI_2_Amber_2PROTEINS2008_2013_08_31.h"

#include "TACCM_CGFGABAbMDrun_Amber_2PROTEINS2008.h"
#include "TACCM_CGAAMDrun_Amber_2PROTEINS2008.h"
#include "TACCM_CGAAMDrun_Amber_2PROTEINS2008_2013_08_31.h"

#include "REMD_functions.h"

#define ON  0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numstep=1000,interval=100;
  int numEX=1,numRE=2;
  int vMode=OFF;

  int equflag=OFF;
  int numstepequ=0,intervalequ=1;

  struct AADataforREMD_Amber FGdata;
  struct CGDataforREMD_PROTEINS2008 CGdata[2];
  struct AACGCommonDataforREMD_A_P2008 Cdata;
  struct TACCMDataforREMD_A_2P2008 Zdata;
  double *KZAA,**KZCG;

  struct potential e;
  struct force f;

  int nc=1;                          
  double T0AA,T0CG[2],T0Z;
  double k_B=1.98723e-3;             
  double UNITT=418.4070;             
  double KTAA,KTCG[2],KBTZ,tau=0.01,tau2,pi;                    

  double fact_b=1.0,fact_a=1.0,fact_t=1.0,fact_NC=1.0,fact_NNC=1.0;

  double dt,dt2,wdt2[3],wdt4[3];

  double **refcrd,ep=0.3,criteria=6.5;
  int NCmode=3,nibnum=3,numnb,num14;
  double EAAm_equ,ECGm_equ[2],EZm_equ;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,**refcrdfilename,*parmfilename,*TACCMfilename;

  char *outputfilenameAAbase,*trjfilenameAAbase,outputfilenameAA[2000],trjfilenameAA[2000];
  char *outputfilenameCG1base,*trjfilenameCG1base,outputfilenameCG1[2000],trjfilenameCG1[2000];
  char *outputfilenameCG2base,*trjfilenameCG2base,outputfilenameCG2[2000],trjfilenameCG2[2000];
  char *trjfileZbase,*trjfileThetaAAbase,*trjfileThetaCG1base,*trjfileThetaCG2base,*logfilename;
  char trjfilenameZ[2000],trjfilenameThetaAA[2000],trjfilenameThetaCG1[2000],trjfilenameThetaCG2[2000],logf[1000];
  FILE *inputfile,**refcrdfile,*parmfile,*TACCMfile,*logfile;

  char *acc_ratio_filename="AccRatio.txt";
  FILE *acc_ratio_file;

  char *progname;
  int opt_idx=1;

  double **acc_ratio;

  int my_rank,num_procs,tag = 0;
  MPI_Status status;            

  MPI_Init(&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  pi=acos(-1.0);

  struct option long_opt[] = {
    {"vMode",1,NULL,'v'},
    {"nums",1,NULL,'s'}, {"int",1,NULL,'i'},
    {"numRE",1,NULL,'n'}, {"numEX",1,NULL,'e'},
    {"tau",1,NULL,'a'}, {"dt",1,NULL,'x'},
    {"TAA",1,NULL,'T'}, {"TCG1",1,NULL,'1'}, {"TCG2",1,NULL,'2'},  {"TZ",1,NULL,'B'},
    {"mZ",1,NULL,'m'}, {"massX",1,NULL,'X'},
    {"ep",1,NULL,'p'}, {"cutoff",1,NULL,'c'},
    {"AR",1,NULL,'A'},
    {"equ",1,NULL,'E'}, {"intequ",1,NULL,'I'},
    {"h",0,NULL,'h'},
    {"fb",1,NULL,'{'},
    {"fa",1,NULL,'}'},
    {"ft",1,NULL,'@'},
    {"fc",1,NULL,'!'},
    {"fn",1,NULL,'>'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hv:n:s:e:a:i:x:T:1:2:B:m:X:p:c:A:E:I:{:}:@:!:>:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 's':
      numstep=atoi(optarg); break;
    case 'i':
      interval=atoi(optarg);  break;
    case 'e':
      numEX=atoi(optarg); break;
    case 'n':
      numRE=atoi(optarg); break;
    case 'T':
      T0AA=atof(optarg);  break;
    case '1':
      T0CG[0]=atof(optarg);  break;
    case '2':
      T0CG[1]=atof(optarg);  break;
    case 'B':
      T0Z=atof(optarg);break;
    case 'v':
      vMode=ON;        break;
    case 'a':
      tau=atof(optarg); break;
    case 'x':
      dt=atof(optarg);   break;
    case 'm':
      Zdata.massZ=atof(optarg);  break;
    case 'p':
      ep=atof(optarg); break;
    case 'c':
      criteria=atof(optarg); break;
    case 'A':
      acc_ratio_filename=optarg;  break;
    case 'E':
      equflag=ON;
      numstepequ=atoi(optarg);  break;
    case 'I':
      intervalequ=atoi(optarg);  break;
    case '{':
      fact_b=atof(optarg);  break;
    case '}':
      fact_a=atof(optarg);  break;
    case '@':
      fact_t=atof(optarg);  break;
    case '!':
      fact_NC=atof(optarg);  break;
    case '>':
      fact_NNC=atof(optarg);  break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);   exit(1);
    }
  }

  progname=*argv;  argc-=optind;  argv+=optind;

  refcrdfilename=(char **)gcemalloc(sizeof(char *)*2);
  for (i=0;i<2;++i) refcrdfilename[i]=(char *)gcemalloc(sizeof(char)*2000);

  if (argc < 16) {
    USAGE(progname);
    exit(1);
  }
  inputfilename        = *argv;
  refcrdfilename[0]    = *++argv;			   
  refcrdfilename[1]    = *++argv;			   
  TACCMfilename        = *++argv;			 
  parmfilename         = *++argv;
  outputfilenameAAbase = *++argv;
  trjfilenameAAbase    = *++argv;
  outputfilenameCG1base= *++argv;
  trjfilenameCG1base   = *++argv;
  outputfilenameCG2base = *++argv;
  trjfilenameCG2base   = *++argv;
  trjfileZbase         = *++argv;
  trjfileThetaAAbase   = *++argv;
  trjfileThetaCG1base  = *++argv;
  trjfileThetaCG2base  = *++argv;
  logfilename          = *++argv;

  index_replicas=(int *)gcemalloc(sizeof(int)*numRE);
  index_parameters=(int *)gcemalloc(sizeof(int)*numRE);
  REMD_ini_purmutation_funcs(numRE);

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 

  Cdata.numatom=AP.NATOM;
  j=0;  for (i=0;i<Cdata.numatom;++i)  if (strncmp(AP.IGRAPH[i],"H",1)==0)  ++j;
  Cdata.numheavyatom=Cdata.numatom-j;  Cdata.numres=AP.NRES;

  TACCMfile=efopen(TACCMfilename,"r");
  fscanf(TACCMfile,"%d",&Zdata.numZ);
  Zdata.pairs=(int **)gcemalloc(sizeof(int *)*Zdata.numZ);
  for (i=0;i<Zdata.numZ;++i) Zdata.pairs[i]=(int *)gcemalloc(sizeof(int)*5);
  for (i=0;i<Zdata.numZ;++i) {
    for (j=0;j<4;++j) fscanf(TACCMfile,"%d",&(Zdata.pairs[i][j]));
    fscanf(TACCMfile,"%d",&(Zdata.pairs[i][j]));
  }
  fclose(TACCMfile);

  Cdata.mass=(double *)gcemalloc(sizeof(double)*Cdata.numatom);
  for (i=0;i<Cdata.numatom;++i) Cdata.mass[i]=AP.AMASS[i];
  refcrd=(double **)gcemalloc(sizeof(double *)*2);
  for (i=0;i<2;++i) refcrd[i]=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);

  refcrdfile=(FILE **)gcemalloc(sizeof(FILE *)*2);

  for (i=0;i<2;++i) {
    refcrdfile[i]=efopen(refcrdfilename[i],"r");
    getline(&line,&len,refcrdfile[i]);
    fscanf(refcrdfile[i],"%d",&d);
    for (j=0;j<Cdata.numatom;++j) for (k=0;k<3;++k) fscanf(refcrdfile[i],"%lf",&refcrd[i][j*3+k]);
    fclose(refcrdfile[i]);
  }

  FGdata.crd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  FGdata.vel=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  CGdata[0].crd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  CGdata[1].crd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  CGdata[0].vel=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  CGdata[1].vel=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);

  Zdata.Z=(double *)gcemalloc(sizeof(double)*Zdata.numZ);
  Zdata.velZ=(double *)gcemalloc(sizeof(double)*Zdata.numZ);

  KZAA=(double *)gcemalloc(sizeof(double)*numRE);
  KZCG=(double **)gcemalloc(sizeof(double *)*2);
  for (i=0;i<2;++i) KZCG[i]=(double *)gcemalloc(sizeof(double)*numRE);
  inputfile=efopen(inputfilename,"r");
  CGAAREMDreadInputs_Amber_PROTEINS2008_Amber_hybrid_1FG2CG(inputfile,Cdata.numatom,numRE,my_rank,
							    FGdata.crd,FGdata.vel,
							    CGdata[0].crd,CGdata[0].vel,
							    CGdata[1].crd,CGdata[1].vel,
							    KZAA,KZCG[0],KZCG[1]);
  fclose(inputfile);
  Zdata.KZAA=KZAA[my_rank];
  for (i=0;i<2;++i) Zdata.KZCG[i]=KZCG[i][my_rank];

  if ( vMode==OFF ) {
    MD_Generate_inivelo(FGdata.vel,Cdata.mass,Cdata.numatom,k_B*T0AA*UNITT);
    MD_Generate_inivelo(CGdata[0].vel,Cdata.mass,Cdata.numatom,k_B*T0CG[0]*UNITT);
    MD_Generate_inivelo(CGdata[1].vel,Cdata.mass,Cdata.numatom,k_B*T0CG[1]*UNITT);
    TACCM_MD_Generate_inivelo(Zdata.velZ,Zdata.massZ,Zdata.numZ,k_B*T0Z*UNITT);

    for (j=0;j<Cdata.numatom;++j) {
      if (strncmp(AP.IGRAPH[j],"H",1)==0) {
	for (k=0;k<3;++k) {
	  CGdata[0].vel[j*3+k]=0.0;
	  CGdata[1].vel[j*3+k]=0.0;
	}
      }
    }

    FGdata.zeta=0.0;      FGdata.V_zeta=0.0;
    CGdata[0].zeta=0.0;      CGdata[1].V_zeta=0.0;
    CGdata[0].zeta=0.0;      CGdata[1].V_zeta=0.0;
    Zdata.zetaZ=0.0;      Zdata.V_zetaZ=0.0;
  }  
  TACCM_CTheta(FGdata.crd,Cdata.numatom,Zdata.Z,Zdata.numZ,Zdata.pairs,pi);

  tau=tau/2.0/pi;  tau2=tau*tau;  
  KBTZ=k_B*T0Z;  Zdata.NfKTZ=(Zdata.numZ+1)*KBTZ*UNITT;  Zdata.QZ=tau2*KBTZ*UNITT*Zdata.numZ;
  KTAA=k_B*T0AA;  FGdata.NfKT=(3.0*Cdata.numatom+1)*KTAA*UNITT;  FGdata.Q=tau2*KTAA*UNITT*(3.0*Cdata.numatom);
  KTCG[0]=k_B*T0CG[0];  CGdata[0].NfKT=(3.0*Cdata.numheavyatom+1)*KTCG[0]*UNITT;  
  CGdata[0].Q=tau2*KTCG[0]*UNITT*(3.0*Cdata.numheavyatom);
  KTCG[1]=k_B*T0CG[1];  CGdata[1].NfKT=(3.0*Cdata.numheavyatom+1)*KTCG[1]*UNITT;  
  CGdata[1].Q=tau2*KTCG[1]*UNITT*(3.0*Cdata.numheavyatom);

  ffL_set_calcffandforce(&(FGdata.e),&(FGdata.f));

  ffL_set_calcffandforce(&e,&f);
  ffL_set_non_bonding_index_1(&numnb,&num14);
  e.parm.numnb=numnb;
  e.parm.num14=num14;
  e.parm.indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  e.parm.index14=(int *)gcemalloc(sizeof(int)*num14*2);
  ffL_set_non_bonding_index_2(e.parm.indexnb,e.parm.index14);

  GOLMAA_PROTEINS2008_ff_set_calcff_b(&(CGdata[0].e),refcrd[0],Cdata.numatom,Cdata.numres,
				      /*FGdata.*/e.parm.indexnb,/*FGdata.*/e.parm.numnb,ep,nibnum,criteria);

  GOLMAA_PROTEINS2008_ff_set_calcff_b(&(CGdata[1].e),refcrd[1],Cdata.numatom,Cdata.numres,
				      /*FGdata.*/e.parm.indexnb,/*FGdata.*/e.parm.numnb,ep,nibnum,criteria);


  MD_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);

  sprintf(outputfilenameAA,"%s_%d",outputfilenameAAbase,my_rank+1);
  sprintf(trjfilenameAA,"%s_%d",trjfilenameAAbase,my_rank+1);

  sprintf(outputfilenameCG1,"%s_%d",outputfilenameCG1base,my_rank+1);
  sprintf(trjfilenameCG1,"%s_%d",trjfilenameCG1base,my_rank+1);

  sprintf(outputfilenameCG2,"%s_%d",outputfilenameCG2base,my_rank+1);
  sprintf(trjfilenameCG2,"%s_%d",trjfilenameCG2base,my_rank+1);

  sprintf(trjfilenameZ,"%s_%d",trjfileZbase,my_rank+1);
  sprintf(trjfilenameThetaAA,"%s_%d",trjfileThetaAAbase,my_rank+1);
  sprintf(trjfilenameThetaCG1,"%s_%d",trjfileThetaCG1base,my_rank+1);
  sprintf(trjfilenameThetaCG2,"%s_%d",trjfileThetaCG2base,my_rank+1);

  sprintf(logf,"%s_%d_ex.log",logfilename,my_rank+1);

  myncL_create_def_MCD(trjfilenameAA,Cdata.numatom,&(FGdata.nc_id_MCD));
  FGdata.outputfile=efopen(outputfilenameAA,"w");

  myncL_create_def_MCD(trjfilenameCG1,Cdata.numatom,&(CGdata[0].nc_id_MCD));
  CGdata[0].outputfile=efopen(outputfilenameCG1,"w");

  myncL_create_def_MCD(trjfilenameCG2,Cdata.numatom,&(CGdata[1].nc_id_MCD));
  CGdata[1].outputfile=efopen(outputfilenameCG2,"w");

  Zdata.trjfileZ=efopen(trjfilenameZ,"w");
  Zdata.trjfilThetaAA=efopen(trjfilenameThetaAA,"w");
  Zdata.trjfilThetaCG1=efopen(trjfilenameThetaCG1,"w");
  Zdata.trjfilThetaCG2=efopen(trjfilenameThetaCG2,"w");

  MPI_Barrier(MPI_COMM_WORLD);

  if ( num_procs != numRE ) {    printf("condition error\n");    exit(1);  }

  logfile=efopen(logf,"w");

  if (equflag==ON) {
    runTACCM_2CG1FG_MD_NHC_MP1998_Amber_PROTEINS2008_2013_08_31(// AA /////////////////////////////////////////////////////////
								FGdata.crd,FGdata.vel,
								&(FGdata.zeta),&(FGdata.V_zeta),FGdata.Q,
								FGdata.e,FGdata.f,FGdata.T,FGdata.NfKT,
								FGdata.avePE,FGdata.aveKE,FGdata.aveT,
								FGdata.varPE,FGdata.varKE,FGdata.varT,
								FGdata.nc_id_MCD,FGdata.outputfile,
						     // CG1 ////////////////////////////////////////////////////////
								CGdata[0].crd,CGdata[0].vel,
								&(CGdata[0].zeta),&(CGdata[0].V_zeta),CGdata[0].Q,
								CGdata[0].e,
								CGdata[0].T,CGdata[0].NfKT,
								CGdata[0].avePE,CGdata[0].aveKE,CGdata[0].aveT,
								CGdata[0].varPE,CGdata[0].varKE,CGdata[0].varT,
								CGdata[0].nc_id_MCD,CGdata[0].outputfile,
						     // CG1 ////////////////////////////////////////////////////////
								CGdata[1].crd,CGdata[1].vel,
								&(CGdata[1].zeta),&(CGdata[1].V_zeta),CGdata[1].Q,
								CGdata[1].e,
								CGdata[1].T,CGdata[1].NfKT,
								CGdata[1].avePE,CGdata[1].aveKE,CGdata[1].aveT,
								CGdata[1].varPE,CGdata[1].varKE,CGdata[1].varT,
								CGdata[1].nc_id_MCD,CGdata[1].outputfile,
								// Z  /////////////////////////////////////////////////////////
								Zdata.Z,Zdata.velZ,Zdata.massZ,
								&(Zdata.zetaZ),&(Zdata.V_zetaZ),
								Zdata.QZ,Zdata.NfKTZ,Zdata.T,
								Zdata.numZ,Zdata.KZAA,Zdata.KZCG[0],Zdata.KZCG[1],Zdata.pairs,
								Zdata.avePEZ,Zdata.aveKEZ,Zdata.aveTZ,
								Zdata.varPEZ,Zdata.varKEZ,Zdata.varTZ, 
								Zdata.trjfileZ,Zdata.trjfilThetaAA,
								Zdata.trjfilThetaCG1,Zdata.trjfilThetaCG2,
								// CM  ////////////////////////////////////////////////////////
								Cdata.mass,(Cdata.numatom),(Cdata.numheavyatom),
								numstepequ,intervalequ,&l,
								dt,dt2,wdt2,wdt4,nc,UNITT,k_B,pi,
								&EAAm_equ,&ECGm_equ[0],&ECGm_equ[1],&EZm_equ,
								fact_b,fact_a,fact_t,fact_NC,fact_NNC);

    FGdata.zeta=0.0; FGdata.V_zeta=0.0;
    CGdata[0].zeta=0.0; CGdata[0].V_zeta=0.0;
    CGdata[1].zeta=0.0; CGdata[1].V_zeta=0.0;
    Zdata.zetaZ=0.0; Zdata.V_zetaZ=0.0;
  }

  acc_ratio=MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_Amber_2PROTEINS2008_2013_08_31(my_rank,num_procs,tag,&status,
										numRE,numEX,KZAA,KZCG,
										FGdata,CGdata,Zdata,Cdata,
										T0AA,T0CG,T0Z,
										numstep,interval,
										dt,dt2,wdt2,wdt4,nc,
										UNITT,k_B,tau,pi,
										logfile,
										fact_b,fact_a,fact_t,fact_NC,fact_NNC);
  
  fclose(logfile);
  fclose(FGdata.outputfile);  
  fclose(CGdata[0].outputfile);  
  fclose(CGdata[1].outputfile);  
  fclose(Zdata.trjfileZ);  
  fclose(Zdata.trjfilThetaAA); 
  fclose(Zdata.trjfilThetaCG1);
  fclose(Zdata.trjfilThetaCG2);

  if ( my_rank==0 ) {
    acc_ratio_file=efopen(acc_ratio_filename,"w");
    for (j=0;j<numRE;++j) {
      if (j!=numRE-1) fprintf(acc_ratio_file,"%d-%d %8.4lf\n",j,j+1,acc_ratio[j][j+1]);
      else fprintf(acc_ratio_file,"%d-0 %8.4lf\n",j,acc_ratio[j][0]);
    }
    fclose(acc_ratio_file);
  }
  
  MPI_Finalize();

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename refcrdfilename TACCMfilename parmfilename outputfilenameAAbase trjfilenameAAbase outputfilenameCGbase trjfilenameCGbase trjfileZbase trjfileThetaZ\n",progname);
}
