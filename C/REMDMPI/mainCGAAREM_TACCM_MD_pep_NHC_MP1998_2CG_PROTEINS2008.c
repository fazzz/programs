
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "REMDCGAA_TACCM_MPI_2_Amber_PROTEINS2008.h"
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

  struct CGDataforREMD_PROTEINS2008 CG1data,CG2data;
  struct AACGCommonDataforREMD_A_P2008 Cdata;
  struct TACCMDataforREMD_A_P2008 Zdata;
  double *KZCG1,*KZCG2;

  struct potential e;
  struct force f;

  int nc=1;                          
  double T0CG1,T0CG2,T0Z;
  double k_B=1.98723e-3;             
  double UNITT=418.4070;             
  double KTCG1,KTCG2,KBTZ,tau=0.01,tau2,pi;                    

  double dt,dt2,wdt2[3],wdt4[3];

  double *refcrd1,*refcrd2,ep=0.3,criteria=6.5;
  int NCmode=3,nibnum=3,numnb,num14;
  double ECG1m_equ,ECG2m_equ,EZm_equ;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*refcrdfilename1,*refcrdfilename2,*parmfilename,*TACCMfilename;

  char *outputfilenameCG1base,*trjfilenameCG1base,outputfilenameCG1[2000],trjfilenameCG1[2000];
  char *outputfilenameCG2base,*trjfilenameCG2base,outputfilenameCG2[2000],trjfilenameCG2[2000];
  char *trjfileZbase,*trjfileThetaCG1base,*trjfileThetaCG2base,*logfilename;
  char trjfilenameZ[2000],trjfilenameThetaCG1[2000],trjfilenameThetaCG2[2000],logf[1000];
  FILE *inputfile,*refcrdfile1,*refcrdfile2,*parmfile,*TACCMfile,*logfile;

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
    {"TCG1",1,NULL,'T'}, {"TCG2",1,NULL,'t'},  {"TZ",1,NULL,'B'},
    {"mZ",1,NULL,'m'}, {"massX",1,NULL,'X'},
    {"ep",1,NULL,'p'}, {"cutoff",1,NULL,'c'},
    {"AR",1,NULL,'A'},
    {"equ",1,NULL,'E'}, {"intequ",1,NULL,'I'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hv:n:s:e:a:i:x:T:t:B:m:X:p:c:A:E:I:",long_opt,&opt_idx))!=-1) {
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
      T0CG1=atof(optarg);  break;
    case 't':
      T0CG2=atof(optarg);  break;
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
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);   exit(1);
    }
  }

  progname=*argv;  argc-=optind;  argv+=optind;

  if (argc < 13) {
    USAGE(progname);
    exit(1);
  }
  inputfilename        = *argv;
  refcrdfilename1      = *++argv;			   
  refcrdfilename2      = *++argv;			   
  TACCMfilename        = *++argv;			 
  parmfilename         = *++argv;
  outputfilenameCG1base = *++argv;
  trjfilenameCG1base    = *++argv;
  outputfilenameCG2base = *++argv;
  trjfilenameCG2base    = *++argv;
  trjfileZbase         = *++argv;
  trjfileThetaCG1base   = *++argv;
  trjfileThetaCG2base   = *++argv;
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
  refcrd1=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  refcrd2=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);

  refcrdfile1=efopen(refcrdfilename1,"r");
  getline(&line,&len,refcrdfile1);
  fscanf(refcrdfile1,"%d",&d);
  for (i=0;i<Cdata.numatom;++i) for (j=0;j<3;++j) fscanf(refcrdfile1,"%lf",&refcrd1[i*3+j]);
  fclose(refcrdfile1);

  refcrdfile2=efopen(refcrdfilename2,"r");
  getline(&line,&len,refcrdfile2);
  fscanf(refcrdfile2,"%d",&d);
  for (i=0;i<Cdata.numatom;++i) for (j=0;j<3;++j) fscanf(refcrdfile2,"%lf",&refcrd2[i*3+j]);
  fclose(refcrdfile2);

  CG1data.crd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  CG1data.vel=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  CG2data.crd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  CG2data.vel=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);

  Zdata.Z=(double *)gcemalloc(sizeof(double)*Zdata.numZ);
  Zdata.velZ=(double *)gcemalloc(sizeof(double)*Zdata.numZ);

  KZCG1=(double *)gcemalloc(sizeof(double)*numRE);
  KZCG2=(double *)gcemalloc(sizeof(double)*numRE);
  inputfile=efopen(inputfilename,"r");
  CGAAREMDreadInputs_test(inputfile,Cdata.numatom,numRE,my_rank,
			  CG1data.crd,CG1data.vel,CG2data.crd,CG2data.vel,
			  KZCG1,KZCG2);
  fclose(inputfile);
  Zdata.KZAA=KZCG1[my_rank];
  Zdata.KZCG=KZCG2[my_rank];

  if ( vMode==OFF ) {
    MD_Generate_inivelo(CG1data.vel,Cdata.mass,Cdata.numatom,k_B*T0CG1*UNITT);
    MD_Generate_inivelo(CG2data.vel,Cdata.mass,Cdata.numatom,k_B*T0CG2*UNITT);
    TACCM_MD_Generate_inivelo(Zdata.velZ,Zdata.massZ,Zdata.numZ,k_B*T0Z*UNITT);

    for (j=0;j<Cdata.numatom;++j) 
      if (strncmp(AP.IGRAPH[j],"H",1)==0)  
	for (k=0;k<3;++k) 
	  CG2data.vel[j*3+k]=0.0;

    CG1data.zeta=0.0;      CG1data.V_zeta=0.0;
    CG2data.zeta=0.0;      CG2data.V_zeta=0.0;
    Zdata.zetaZ=0.0;      Zdata.V_zetaZ=0.0;
  }  
  TACCM_CTheta(CG1data.crd,Cdata.numatom,Zdata.Z,Zdata.numZ,Zdata.pairs,pi);

  tau=tau/2.0/pi;  tau2=tau*tau;  
  KBTZ=k_B*T0Z;  Zdata.NfKTZ=(Zdata.numZ+1)*KBTZ*UNITT;  Zdata.QZ=tau2*KBTZ*UNITT*Zdata.numZ;
  KTCG1=k_B*T0CG1;  CG1data.NfKT=(3.0*Cdata.numatom+1)*KTCG1*UNITT;  CG1data.Q=tau2*KTCG1*UNITT*(3.0*Cdata.numatom);
  KTCG2=k_B*T0CG2;  CG2data.NfKT=(3.0*Cdata.numheavyatom+1)*KTCG2*UNITT;  
  CG2data.Q=tau2*KTCG2*UNITT*(3.0*Cdata.numheavyatom);

  //  ffL_set_calcffandforce(&(CG1data.e),&(CG1data.f));

  ffL_set_calcffandforce(&e,&f);
  ffL_set_non_bonding_index_1(&numnb,&num14);
  e.parm.numnb=numnb;
  e.parm.num14=num14;
  e.parm.indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  e.parm.index14=(int *)gcemalloc(sizeof(int)*num14*2);
  ffL_set_non_bonding_index_2(e.parm.indexnb,e.parm.index14);

  GOLMAA_PROTEINS2008_ff_set_calcff_b(&(CG1data.e),refcrd1,Cdata.numatom,Cdata.numres,
  				      /*AAdata.*/e.parm.indexnb,/*AAdata.*/e.parm.numnb,ep,nibnum,criteria);

  ffL_set_calcffandforce(&e,&f);
  ffL_set_non_bonding_index_1(&numnb,&num14);
  e.parm.numnb=numnb;
  e.parm.num14=num14;
  e.parm.indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  e.parm.index14=(int *)gcemalloc(sizeof(int)*num14*2);
  ffL_set_non_bonding_index_2(e.parm.indexnb,e.parm.index14);

  GOLMAA_PROTEINS2008_ff_set_calcff_b(&(CG2data.e),refcrd2,Cdata.numatom,Cdata.numres,
  				      /*AAdata.*/e.parm.indexnb,/*AAdata.*/e.parm.numnb,ep,nibnum,criteria);

  MD_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);

  sprintf(outputfilenameCG1,"%s_%d",outputfilenameCG1base,my_rank+1);
  sprintf(trjfilenameCG1,"%s_%d",trjfilenameCG1base,my_rank+1);

  sprintf(outputfilenameCG2,"%s_%d",outputfilenameCG2base,my_rank+1);
  sprintf(trjfilenameCG2,"%s_%d",trjfilenameCG2base,my_rank+1);

  sprintf(trjfilenameZ,"%s_%d",trjfileZbase,my_rank+1);
  sprintf(trjfilenameThetaCG1,"%s_%d",trjfileThetaCG1base,my_rank+1);
  sprintf(trjfilenameThetaCG2,"%s_%d",trjfileThetaCG2base,my_rank+1);

  sprintf(logf,"%s_%d_ex.log",logfilename,my_rank+1);

  myncL_create_def_MCD(trjfilenameCG1,Cdata.numatom,&(CG1data.nc_id_MCD));
  CG1data.outputfile=efopen(outputfilenameCG1,"w");

  myncL_create_def_MCD(trjfilenameCG2,Cdata.numatom,&(CG2data.nc_id_MCD));
  CG2data.outputfile=efopen(outputfilenameCG2,"w");

  Zdata.trjfileZ=efopen(trjfilenameZ,"w");
  Zdata.trjfilThetaAA=efopen(trjfilenameThetaCG1,"w");
  Zdata.trjfilThetaCG=efopen(trjfilenameThetaCG2,"w");

  MPI_Barrier(MPI_COMM_WORLD);

  if ( num_procs != numRE ) {    printf("condition error\n");    exit(1);  }

  logfile=efopen(logf,"w");

  if (equflag==ON) {
    runTACCM_CGAA_MD_NHC_MP1998_2PROTEINS2008(// CG1 /////////////////////////////////////////////////////////
					      CG1data.crd,CG1data.vel,
					      &(CG1data.zeta),&(CG1data.V_zeta),CG1data.Q,
					      CG1data.e,
					      CG1data.T,CG1data.NfKT,
					      CG1data.avePE,CG1data.aveKE,CG1data.aveT,
					      CG1data.varPE,CG1data.varKE,CG1data.varT,
					      CG1data.nc_id_MCD,CG1data.outputfile,
					      // CG2 /////////////////////////////////////////////////////////
					      CG2data.crd,CG2data.vel,
					      &(CG2data.zeta),&(CG2data.V_zeta),CG2data.Q,
					      CG2data.e,
					      CG2data.T,CG2data.NfKT,
					      CG2data.avePE,CG2data.aveKE,CG2data.aveT,
					      CG2data.varPE,CG2data.varKE,CG2data.varT,
					      CG2data.nc_id_MCD,CG2data.outputfile,
					      // Z  /////////////////////////////////////////////////////////
					      Zdata.Z,Zdata.velZ,Zdata.massZ,
					      &(Zdata.zetaZ),&(Zdata.V_zetaZ),
					      Zdata.QZ,Zdata.NfKTZ,Zdata.T,
					      Zdata.numZ,Zdata.KZAA,Zdata.KZCG,Zdata.pairs,
					      Zdata.avePEZ,Zdata.aveKEZ,Zdata.aveTZ,
					      Zdata.varPEZ,Zdata.varKEZ,Zdata.varTZ, 
					      Zdata.trjfileZ,Zdata.trjfilThetaAA,Zdata.trjfilThetaCG,
					      // CM  /////////////////////////////////////////////////////////
					      Cdata.mass,(Cdata.numatom),(Cdata.numheavyatom),
					      numstepequ,intervalequ,&l,
					      dt,dt2,wdt2,wdt4,nc,UNITT,k_B,pi,
					      &ECG1m_equ,&ECG2m_equ,&EZm_equ);

    CG1data.zeta=0.0; CG1data.V_zeta=0.0;
    CG2data.zeta=0.0; CG2data.V_zeta=0.0;
    Zdata.zetaZ=0.0; Zdata.V_zetaZ=0.0;

  }

  acc_ratio=MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_2PROTEINS2008(my_rank,num_procs,tag,&status,
							       numRE,numEX,KZCG1,KZCG2,
							       CG1data,CG2data,Zdata,Cdata,
							       T0CG1,T0CG2,T0Z,
							       numstep,interval,
							       dt,dt2,wdt2,wdt4,nc,
							       UNITT,k_B,tau,pi,
							       logfile);
  
  fclose(logfile);
  fclose(CG1data.outputfile);  
  fclose(CG2data.outputfile);  
  fclose(Zdata.trjfileZ);  
  fclose(Zdata.trjfilThetaAA); 
  fclose(Zdata.trjfilThetaCG);

  if ( my_rank==0 ) {
    acc_ratio_file=efopen(acc_ratio_filename,"w");
    for (i=0;i<numRE;++i) {
      //      for (j=i+1;j<numRE;++j) {
      if (i!=numRE-1) fprintf(acc_ratio_file,"%d-%d %8.4lf\n",i,i+1,acc_ratio[i][i+1]);
      else fprintf(acc_ratio_file,"%d-0 %8.4lf\n",i,acc_ratio[i][0]);
	//      }
    }
    fclose(acc_ratio_file);
  }
  
  MPI_Finalize();

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename refcrdfilename TACCMfilename parmfilename outputfilenameCG1base trjfilenameCG1base outputfilenameCG2base trjfilenameCG2base trjfileZbase trjfileThetaZ\n",progname);
}
