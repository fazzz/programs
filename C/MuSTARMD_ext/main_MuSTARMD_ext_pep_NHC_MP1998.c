
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "MuSTARMD_ext_MPI.h"
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

  int ND=2;

  struct AAData_MuSTARMD *Alphadata;
  struct TACCMData_MuSTARMD Zdata;
  struct AACGCommonData_MuSTARMD  Cdata;

  struct AmberParmL *ap_Alpha;

  double **KZAlpha;

  int nc=1;                          
  double *T0Alpha,T0Z;
  double k_B=1.98723e-3;             
  double UNITT=418.4070;             
  double *KTAlpha,KBTZ,tau=0.01,tau2,pi;

  double dt,dt2,wdt2[3],wdt4[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,**parmfilenameAlpha,*TACCMfilename;

  char **outputfilenameAlphabase,**trjfilenameAlphabase;
  char **outputfilenameAlpha,**trjfilenameAlpha;
  char *trjfileZbase,**trjfileThetaAlphabase,*logfilename;
  char *trjfilenameZ,**trjfilenameThetaAlpha,*logf;
  FILE *inputfile,**parmfileAlpha,*TACCMfile,*logfile;

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
    {"numD",1,NULL,'D'},
    {"nums",1,NULL,'s'}, {"int",1,NULL,'i'}, {"numEX",1,NULL,'e'},
    {"tau",1,NULL,'a'}, {"dt",1,NULL,'x'},
    {"mZ",1,NULL,'m'}, {"AR",1,NULL,'A'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hv:D:s:i:e:a:x:m:A:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'v':
      vMode=ON;        break;
    case 'D':
      numD=atoi(optarg); break;
    case 's':
      numstep=atoi(optarg); break;
    case 'i':
      interval=atoi(optarg);  break;
    case 'e':
      numEX=atoi(optarg); break;
    case 'a':
      tau=atof(optarg); break;
    case 'x':
      dt=atof(optarg);   break;
    case 'm':
      Zdata.massZ=atof(optarg);  break;
    case 'p':
      ep=atof(optarg); break;
    case 'A':
      acc_ratio_filename=optarg;  break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);   exit(1);
    }
  }

  parmfilenameAlpha=(char **)gcemalloc(sizeof(char *)*ND);
  outputfilenameAlphabase=(char **)gcemalloc(sizeof(char *)*ND);
  trjfilenameAlphabase=(char **)gcemalloc(sizeof(char *)*ND);
  trjfileThetaAlphabase=(char **)gcemalloc(sizeof(char *)*ND);
  T0Alpha=(double **)gcemalloc(sizeof(double *)*ND);

  progname=*argv;  argc-=optind;  argv+=optind;

  if (argc < 5+5*ND) {
    USAGE(progname);
    exit(1);
  }
  inputfilename        = *argv;
  TACCMfilename        = *++argv;			 
  for (i=0;i<ND;++i) parmfilenameAlpha[i]       = *++argv;
  for (i=0;i<ND;++i) outputfilenameAlphabase[i] = *++argv;
  for (i=0;i<ND;++i) trjfilenameAlphabase[i]    = *++argv;
  for (i=0;i<ND;++i) T0Alpha[i]           = atof(*++argv);
  T0Z                  = atof(*++argv);
  trjfileZbase         = *++argv;
  for (i=0;i<ND;++i) trjfileThetaAlphabase[i]   = *++argv;
  logfilename          = *++argv;

  index_replicas=(int *)gcemalloc(sizeof(int)*numRE);
  index_parameters=(int *)gcemalloc(sizeof(int)*numRE);
  REMD_ini_purmutation_funcs(numRE);

  for (i=0i<ND;++i) {
    parmfileAlpha[i]=efopen(parmfilenameAlpha[i],"r");
    readParmtopLb(parmfileAlpha[i],&ap_Alpha[i]);
    fclose(parmfileAlpha[i]); 
  }

  for (i=1;i<ND;++i) {
    if ( ap_Alpha[0].NATOM != ap_Alpha[1].NATOM ) {
      printf("error: about num of atoms\n");
    }
  }

  Cdata.numatom=ap_Alpha[0].NATOM;

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

  for (i=0;i<ND;++i) {
    Alpa[i].crd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3*2);
    Alpa[i].vel=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3*2);
  }
  Zdata.Z=(double *)gcemalloc(sizeof(double)*Zdata.numZ);
  Zdata.velZ=(double *)gcemalloc(sizeof(double)*Zdata.numZ);

  KZAlpa=(double **)gcemalloc(sizeof(double *));
  for (i=0;i<ND;++i) KZAlpa[i]=(double *)gcemalloc(sizeof(double)*numRE);

  inputfile=efopen(inputfilename,"r");
  MuSTARMD_ext_readInputs(inputfile,ND,
			  Cdata.numatom,numRE,my_rank,
			  Alphadata,KZAlpha);
  fclose(inputfile);

  for (i=0;i<ND;++i) {
    Zdata.KZAlpa=KZAlpa[my_rank];
  }

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

  ffLc_set_calcffandforce(&(AAdata.e),&(AAdata.f),ap_AA);
  ffLc_set_calcffandforce(&(CGdata.e),&(CGdata.f),ap_CG);

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

  logfile=efopen(logf,"w");
  acc_ratio=MPI_MuSTARMD_ext_pep_NHC_MP1998(my_rank,num_procs,tag,&status,
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
  printf("%s [-h]   inputfilename TACCMfilename  parmfilenameAA  parmfilenameCG outputfilenameAAbase trjfilenameAAbase outputfilenameCGbase trjfilenameCGbase trjfileZbase trjfileThetaAAbase trjfileThetaCGbase logfilename\n",progname);
}
