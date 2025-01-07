
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "REMDCGAA_TACCM_MPI_2_testb.h"
#include "REMDCGAA_TACCM_MPI_2_testc.h"
#include "REMD_functions.h"

#include "EXT_SYS.h"

#define ON  0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numstep=1000,interval=100;
  int numEX=1,numRE=2;
  int vMode=OFF;
  int numEXT=1;

  int equflag=OFF;
  int numstepequ=0,intervalequ=1;

  struct AADataforREMD_testb AAdata, CGdata;
  struct AACGCommonDataforREMD_testb Cdata;
  struct TACCMDataforREMD_testb Zdata;

  struct AmberParmL ap_AA,ap_CG;

  double *KZAA,*KZCG;
  double parameterCG=0.2;

  int nc=1;                          
  double T0AA,T0CG,T0Z;
  double k_B=1.98723e-3;             
  double UNITT=418.4070;             
  double KTAA,KTCG,KBTZ,tau=0.01,tau2,pi;                    

  double dt,dt2,wdt2[3],wdt4[3];

  double *refcrd,ep=0.3,criteria=6.5;
  int NCmode=3,nibnum=3,numnb,num14;
  double ECG1m_equ,ECG2m_equ,EZm_equ;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*refcrdfilename,*parmfilenameAA,*parmfilenameCG,*TACCMfilename;

  char *outputfilenameAAbase,*trjfilenameAAbase,outputfilenameAA[2000],trjfilenameAA[2000];
  char *outputfilenameCGbase,*trjfilenameCGbase,outputfilenameCG[2000],trjfilenameCG[2000];
  char *trjfileZbase,*trjfileThetaAAbase,*trjfileThetaCGbase,*logfilename;
  char trjfilenameZ[2000],trjfilenameThetaAA[2000],trjfilenameThetaCG[2000],logf[1000];
  FILE *inputfile,*refcrdfile,*parmfileAA,*parmfileCG,*TACCMfile,*logfile;

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
    {"numEXT",1,NULL,'E'}, {"pCG",1,NULL,'P'},
    {"tau",1,NULL,'a'}, {"dt",1,NULL,'x'},
    {"TAA",1,NULL,'T'}, {"TCG",1,NULL,'t'},  {"TZ",1,NULL,'B'},
    {"mZ",1,NULL,'m'}, {"massX",1,NULL,'X'},
    {"ep",1,NULL,'p'}, {"cutoff",1,NULL,'c'},
    {"AR",1,NULL,'A'},
    {"equ",1,NULL,'F'}, {"intequ",1,NULL,'I'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hv:n:E:s:e:a:i:x:T:t:B:m:X:p:c:E:P:A:F:I:",long_opt,&opt_idx))!=-1) {
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
    case 't':
      T0CG=atof(optarg);  break;
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
    case 'E':
      numEXT=atoi(optarg);  break;
    case 'P':
      parameterCG=atof(optarg);  break;
    case 'A':
      acc_ratio_filename=optarg;  break;
    case 'F':
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
  refcrdfilename       = *++argv;			   
  TACCMfilename        = *++argv;			 
  parmfilenameAA       = *++argv;
  parmfilenameCG       = *++argv;
  outputfilenameAAbase = *++argv;
  trjfilenameAAbase    = *++argv;
  outputfilenameCGbase = *++argv;
  trjfilenameCGbase    = *++argv;
  trjfileZbase         = *++argv;
  trjfileThetaAAbase   = *++argv;
  trjfileThetaCGbase   = *++argv;
  logfilename          = *++argv;

  index_replicas=(int *)gcemalloc(sizeof(int)*numRE);
  index_parameters=(int *)gcemalloc(sizeof(int)*numRE);
  REMD_ini_purmutation_funcs(numRE);

  parmfileAA=efopen(parmfilenameAA,"r");
  readParmtopLb(parmfileAA,&ap_AA);
  fclose(parmfileAA); 

  parmfileCG=efopen(parmfilenameCG,"r");
  readParmtopLb(parmfileCG,&ap_CG);
  fclose(parmfileCG); 

  if ( ap_AA.NATOM != ap_CG.NATOM ) {
    printf("error: about num of atoms\n");
  }

  /******************************************************/
  /* for (i=0;i<ap_AA.NATOM;++i) ap_AA.CHRG[i]=0.0;     */
  /* for (i=0;i<(ap_AA.NTYPES)*(ap_AA.NTYPES+1)/2;++i){ */
  /*   ap_AA.CN1[i]=0.0; ap_AA.CN2[i]=0.0;	        */
  /* }						        */
  /* 						        */
  /* for (i=0;i<ap_CG.NATOM;++i) ap_CG.CHRG[i]=0.0;     */
  /* for (i=0;i<(ap_CG.NTYPES)*(ap_CG.NTYPES+1)/2;++i){ */
  /*   ap_CG.CN1[i]=0.0; ap_CG.CN2[i]=0.0;	        */
  /* }						        */
  /******************************************************/

  Cdata.numatom=ap_AA.NATOM;
  j=0;  for (i=0;i<Cdata.numatom;++i)  if (strncmp(ap_AA.IGRAPH[i],"H",1)==0)  ++j;
  Cdata.numheavyatom=Cdata.numatom-j;  Cdata.numres=ap_AA.NRES;

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
  for (i=0;i<Cdata.numatom;++i) Cdata.mass[i]=ap_AA.AMASS[i];
  refcrd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);

  AAdata.crd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3*2);
  AAdata.vel=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3*2);
  CGdata.crd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3*2);
  CGdata.vel=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3*2);

  Zdata.Z=(double *)gcemalloc(sizeof(double)*Zdata.numZ);
  Zdata.velZ=(double *)gcemalloc(sizeof(double)*Zdata.numZ);

  KZAA=(double *)gcemalloc(sizeof(double)*numRE);
  KZCG=(double *)gcemalloc(sizeof(double)*numRE);
  inputfile=efopen(inputfilename,"r");
  CGAAREMDreadInputs_testc(inputfile,Cdata.numatom/numEXT,numRE,my_rank,
			   AAdata.crd,AAdata.vel,CGdata.crd,CGdata.vel,
			   KZAA,KZCG);
  fclose(inputfile);
  Zdata.KZAA=KZAA[my_rank];
  Zdata.KZCG=KZCG[my_rank];

  /********************************************************************************************/
  /* refcrdfile=efopen(refcrdfilename,"r");						      */
  /* getline(&line,&len,refcrdfile);  fscanf(refcrdfile,"%d",&d);			      */
  /* for (i=0;i<Cdata.numatom;++i) for (j=0;j<3;++j) fscanf(refcrdfile,"%lf",&refcrd[i*3+j]); */
  /* fclose(refcrdfile);								      */
  /********************************************************************************************/

  if ( vMode==OFF ) {
    for (i=0;i<numRE;++i) {
      MD_Generate_inivelo(AAdata.vel,Cdata.mass,Cdata.numatom,k_B*T0AA*UNITT);
      MD_Generate_inivelo(CGdata.vel,Cdata.mass,Cdata.numatom,k_B*T0CG*UNITT);
      TACCM_MD_Generate_inivelo(Zdata.velZ,Zdata.massZ,Zdata.numZ,k_B*T0Z*UNITT);

      for (j=0;j<Cdata.numatom/2;++j) {
	for (k=0;k<3;++k) {
	  CGdata.vel[(j+Cdata.numatom/2)*3+k]=0.0;
	  CGdata.vel[j*3+k]=0.0;
	}
      }

      AAdata.zeta=0.0;      AAdata.V_zeta=0.0;
      CGdata.zeta=0.0;      CGdata.V_zeta=0.0;
      Zdata.zetaZ=0.0;      Zdata.V_zetaZ=0.0;
    }
  }  
  TACCM_CTheta(AAdata.crd,Cdata.numatom,Zdata.Z,Zdata.numZ,Zdata.pairs,pi);

  tau=tau/2.0/pi;  tau2=tau*tau;  
  KBTZ=k_B*T0Z;  Zdata.NfKTZ=(Zdata.numZ+1)*KBTZ*UNITT;  Zdata.QZ=tau2*KBTZ*UNITT*Zdata.numZ;
  KTAA=k_B*T0AA;  AAdata.NfKT=(3.0*Cdata.numatom+1)*KTAA*UNITT;  AAdata.Q=tau2*KTAA*UNITT*(3.0*Cdata.numatom);
  KTCG=k_B*T0CG;  CGdata.NfKT=(3.0*Cdata.numheavyatom+1)*KTCG*UNITT;  
  CGdata.Q=tau2*KTCG*UNITT*(3.0*Cdata.numheavyatom);

  /******************************************************/
  /* for (i=0;i<ap_AA.NATOM;++i) ap_AA.CHRG[i]=0.0;     */
  /* for (i=0;i<(ap_AA.NTYPES)*(ap_AA.NTYPES+1)/2;++i){ */
  /*   ap_AA.CN1[i]=0.0; ap_AA.CN2[i]=0.0;	        */
  /* }						        */
  /* 						        */
  /* for (i=0;i<ap_CG.NATOM;++i) ap_CG.CHRG[i]=0.0;     */
  /* for (i=0;i<(ap_CG.NTYPES)*(ap_CG.NTYPES+1)/2;++i){ */
  /*   ap_CG.CN1[i]=0.0; ap_CG.CN2[i]=0.0;	        */
  /* }						        */
  /******************************************************/

  ffLc_set_calcffandforce(&(AAdata.e),&(AAdata.f),ap_AA);
  ffLc_set_calcffandforce(&(CGdata.e),&(CGdata.f),ap_CG);

  AAdata.e.p_e=(double *)gcemalloc(sizeof(double)*Cdata.numatom);
  AAdata.e.p_LJ=(double *)gcemalloc(sizeof(double)*Cdata.numatom);

  CGdata.e.p_e=(double *)gcemalloc(sizeof(double)*Cdata.numatom);
  CGdata.e.p_LJ=(double *)gcemalloc(sizeof(double)*Cdata.numatom);

  ffLc_calcffandforce(AAdata.crd,Cdata.numatom,&(AAdata.e),&(AAdata.f),ap_AA);
  ffLc_calcffandforce(CGdata.crd,Cdata.numatom,&(CGdata.e),&(CGdata.f),ap_CG);
  //  GOLMAA_PROTEINS2008_ff_set_calcff_b(&(CGdata.e_GOLM),refcrd,Cdata.numatom,Cdata.numres,
  //				      AAdata.e.parm.indexnb,AAdata.e.parm.numnb,ep,nibnum,criteria);

  MDb_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);

  sprintf(outputfilenameAA,"%s_%d",outputfilenameAAbase,my_rank+1);
  sprintf(trjfilenameAA,"%s_%d",trjfilenameAAbase,my_rank+1);

  sprintf(outputfilenameCG,"%s_%d",outputfilenameCGbase,my_rank+1);
  sprintf(trjfilenameCG,"%s_%d",trjfilenameCGbase,my_rank+1);

  sprintf(trjfilenameZ,"%s_%d",trjfileZbase,my_rank+1);
  sprintf(trjfilenameThetaAA,"%s_%d",trjfileThetaAAbase,my_rank+1);
  sprintf(trjfilenameThetaCG,"%s_%d",trjfileThetaCGbase,my_rank+1);

  sprintf(logf,"%s_%d_ex.log",logfilename,my_rank+1);

  myncL_create_def_MCD(trjfilenameAA,Cdata.numatom,&(AAdata.nc_id_MCD));
  AAdata.outputfile=efopen(outputfilenameAA,"w");

  myncL_create_def_MCD(trjfilenameCG,Cdata.numatom,&(CGdata.nc_id_MCD));
  CGdata.outputfile=efopen(outputfilenameCG,"w");

  Zdata.trjfileZ=efopen(trjfilenameZ,"w");
  Zdata.trjfilThetaAA=efopen(trjfilenameThetaAA,"w");
  Zdata.trjfilThetaCG=efopen(trjfilenameThetaCG,"w");

  MPI_Barrier(MPI_COMM_WORLD);

  if ( num_procs != numRE ) {    printf("condition error\n");    exit(1);  }

  for (i=0;i<Cdata.numatom;++i) {
    for (j=0;j<3;++j) {
      AAdata.f.f_b[i*3+j]=0.0;
      AAdata.f.f_a[i*3+j]=0.0;
      AAdata.f.f_d[i*3+j]=0.0;
      AAdata.f.f_e[i*3+j]=0.0;
      AAdata.f.f_LJ[i*3+j]=0.0;
      AAdata.f.f_e_14[i*3+j]=0.0;
      AAdata.f.f_LJ_14[i*3+j]=0.0;

      CGdata.f.f_b[i*3+j]=0.0;
      CGdata.f.f_a[i*3+j]=0.0;
      CGdata.f.f_d[i*3+j]=0.0;
      CGdata.f.f_e[i*3+j]=0.0;
      CGdata.f.f_LJ[i*3+j]=0.0;
      CGdata.f.f_e_14[i*3+j]=0.0;
      CGdata.f.f_LJ_14[i*3+j]=0.0;
    }
  }

  if ( equflag == ON ) {
    runTACCM_CGAA_MDc_NHC_MP1998_test(// AA /////////////////////////////////////////////////////////
				      AAdata.crd,AAdata.vel,
				      &(AAdata.zeta),&(AAdata.V_zeta),AAdata.Q,
				      AAdata.e,AAdata.f,ap_AA,AAdata.T,AAdata.NfKT,
				      AAdata.avePE,AAdata.aveKE,AAdata.aveT,
				      AAdata.varPE,AAdata.varKE,AAdata.varT,
				      AAdata.nc_id_MCD,AAdata.outputfile,
				      // CG /////////////////////////////////////////////////////////
				      CGdata.crd,CGdata.vel,
				      &(CGdata.zeta),&(CGdata.V_zeta),CGdata.Q,
				      CGdata.e,CGdata.f,ap_CG,
				      CGdata.T,CGdata.NfKT,
				      CGdata.avePE,CGdata.aveKE,CGdata.aveT,
				      CGdata.varPE,CGdata.varKE,CGdata.varT,
				      CGdata.nc_id_MCD,CGdata.outputfile,
				      // Z  /////////////////////////////////////////////////////////
				      Zdata.Z,Zdata.velZ,Zdata.massZ,
				      &(Zdata.zetaZ),&(Zdata.V_zetaZ),
				      Zdata.QZ,Zdata.NfKTZ,Zdata.T,
				      Zdata.numZ,Zdata.KZAA,Zdata.KZCG,Zdata.pairs,
				      Zdata.avePEZ,Zdata.aveKEZ,Zdata.aveTZ,
				      Zdata.varPEZ,Zdata.varKEZ,Zdata.varTZ, 
				      Zdata.trjfileZ,Zdata.trjfilThetaAA,Zdata.trjfilThetaCG,
				      // CM  /////////////////////////////////////////////////////////
				      Cdata.mass,(Cdata.numatom),
				      numstepequ,intervalequ,&l,
				      dt,dt2,wdt2,wdt4,nc,UNITT,k_B,pi,
				      &ECG1m_equ,&ECG2m_equ,&EZm_equ);
  }

  AAdata.zeta=0.0; CGdata.V_zeta=0.0;
  Zdata.zetaZ=0.0; Zdata.V_zetaZ=0.0;

  logfile=efopen(logf,"w");
  acc_ratio=MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_testc(my_rank,num_procs,tag,&status,
						       numRE,numEX,KZAA,KZCG,
						       AAdata,CGdata,Zdata,Cdata,
						       ap_AA,ap_CG,
						       T0AA,T0CG,T0Z,
						       numstep,interval,
						       dt,dt2,wdt2,wdt4,nc,
						       UNITT,k_B,tau,pi,
						       parameterCG,logfile);
  
  fclose(logfile);
  fclose(AAdata.outputfile);  
  fclose(CGdata.outputfile);  
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
  printf("%s [-h] inputfilename refcrdfilename TACCMfilename parmfilename outputfilenameAAbase trjfilenameAAbase outputfilenameCGbase trjfilenameCGbase trjfileZbase trjfileThetaZ\n",progname);
}
