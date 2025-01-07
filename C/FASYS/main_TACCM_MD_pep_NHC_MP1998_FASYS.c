#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "REMDCGAA_TACCM_MPI_2_test.h"

#include "EXT_SYS.h"

#define ON  0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d;
  int vMode=OFF;
  int numEXT=1;

  struct AADataforREMD_test AAdata;
  struct AACGCommonDataforREMD_test Cdata;
  struct TACCMDataforREMD_test Zdata;

  double avePE,aveKE,aveT,varPE,varKE,varT;
  double avePEZ,aveKEZ,aveTZ,varPEZ,varKEZ,varTZ;

  double parameterCG=0.2;

  int nc=1;                          
  double T0AA,T0Z;
  double k_B=1.98723e-3;             
  double UNITT=418.4070;             
  double KTAA,KTCG,KBTZ,tau=0.01,tau2,pi;                    

  double dt,dt2,wdt2[3],wdt4[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*parmfilename,*TACCMfilename;

  char *outputfilenameAA,*trjfilenameAA;
  char *trjfilenameZ,*trjfilenameThetaAA;

  FILE *inputfile,*parmfile,*TACCMfile,*outputfile;
  FILE *trjfileZ,*trjfileThetaAA;

  char *progname;
  int opt_idx=1;

  pi=acos(-1.0);

  struct option long_opt[] = {
    {"vMode",1,NULL,'v'},
    {"nums",1,NULL,'s'}, {"int",1,NULL,'i'},
    {"tau",1,NULL,'a'}, {"dt",1,NULL,'x'},
    {"TAA",1,NULL,'T'}, {"TZ",1,NULL,'B'},    
    {"numEXT",1,NULL,'E'},
    {"KZ",1,NULL,'K'},
    {"mZ",1,NULL,'m'}, {"massX",1,NULL,'X'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hv:n:s:a:i:x:T:B:m:E:K:X:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 's':
      Cdata.numstep=atoi(optarg); break;
    case 'i':
      Cdata.interval=atoi(optarg);  break;
    case 'T':
      T0AA=atof(optarg);  break;
    case 'B':
      T0Z=atof(optarg);break;
    case 'v':
      vMode=ON;        break;
    case 'E':
      numEXT=atoi(optarg);  break;
    case 'a':
      tau=atof(optarg); break;
    case 'x':
      dt=atof(optarg);   break;
    case 'K':
      Zdata.KZAA=atof(optarg);  break;
    case 'm':
      Zdata.massZ=atof(optarg);  break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);   exit(1);
    }
  }

  progname=*argv;  argc-=optind;  argv+=optind;

  if (argc < 7) {
    USAGE(progname);
    exit(1);
  }
  inputfilename        = *argv;
  TACCMfilename        = *++argv;			 
  parmfilename         = *++argv;
  outputfilenameAA = *++argv;
  trjfilenameAA    = *++argv;
  trjfilenameZ      = *++argv;
  trjfilenameThetaAA   = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 

  extend_test_system(numEXT);
  for (i=0;i<AP.NATOM;++i) AP.CHRG[i]=0.0;
  for (i=0;i<(AP.NTYPES)*(AP.NTYPES+1)/2;++i){
    AP.CN1[i]=0.0; AP.CN2[i]=0.0;
  }

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

  AAdata.crd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  AAdata.vel=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);

  Zdata.Z=(double *)gcemalloc(sizeof(double)*Zdata.numZ);
  Zdata.velZ=(double *)gcemalloc(sizeof(double)*Zdata.numZ);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);  fscanf(inputfile,"%d",&d);
  for (i=0;i<Cdata.numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&(AAdata.crd[i*3+j]));
  fclose(inputfile);

  for (i=1;i<numEXT;++i)
    for (j=0;j<Cdata.numatom/numEXT*3;++j)
      AAdata.crd[i*Cdata.numatom/numEXT*3+j]=AAdata.crd[j];

  if ( vMode==OFF ) {
    MD_Generate_inivelo(AAdata.vel,Cdata.mass,Cdata.numatom,k_B*T0AA*UNITT);
    TACCM_MD_Generate_inivelo(Zdata.velZ,Zdata.massZ,Zdata.numZ,k_B*T0Z*UNITT);

    AAdata.zeta=0.0;      AAdata.V_zeta=0.0;
    Zdata.zetaZ=0.0;      Zdata.V_zetaZ=0.0;
  }  
  TACCM_CTheta(AAdata.crd,Cdata.numatom,Zdata.Z,Zdata.numZ,Zdata.pairs,pi);

  tau=tau/2.0/pi;  tau2=tau*tau;  
  KBTZ=k_B*T0Z;  Zdata.NfKTZ=(Zdata.numZ+1)*KBTZ*UNITT;  Zdata.QZ=tau2*KBTZ*UNITT*Zdata.numZ;
  KTAA=k_B*T0AA;  AAdata.NfKT=(3.0*Cdata.numatom+1)*KTAA*UNITT;  AAdata.Q=tau2*KTAA*UNITT*(3.0*Cdata.numatom);

  ffL_set_calcffandforce(&(AAdata.e),&(AAdata.f));

  MD_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);

  myncL_create_def_MCD(trjfilenameAA,Cdata.numatom,&(AAdata.nc_id_MCD));
  AAdata.outputfile=efopen(outputfilenameAA,"w");

  Zdata.trjfileZ=efopen(trjfilenameZ,"w");
  Zdata.trjfilThetaAA=efopen(trjfilenameThetaAA,"w");

  runTACCM_MD_NHC_MP1998_Amber_AAFF_test(AAdata.crd,AAdata.vel,Cdata.mass,(Cdata.numatom),
					 &(AAdata.zeta),&(AAdata.V_zeta),AAdata.Q,
					 AAdata.e,AAdata.f,
					 AAdata.T,AAdata.NfKT,Cdata.numstep,Cdata.interval,&l,
					 dt,dt2,wdt2,wdt4,nc,
					 &avePE,&aveKE,&aveT,&varPE,&varKE,&varT,
					 UNITT,k_B,
					 AAdata.nc_id_MCD,AAdata.outputfile,
					 Zdata.Z,Zdata.velZ,Zdata.massZ,
					 &(Zdata.zetaZ),&(Zdata.V_zetaZ),
					 Zdata.T,Zdata.QZ,Zdata.NfKTZ,Zdata.numZ,
					 Zdata.KZAA,Zdata.pairs,pi,
					 &avePEZ,&aveKEZ,&aveTZ,&varPEZ,&varKEZ,&varTZ, 
					 Zdata.trjfileZ,Zdata.trjfilThetaAA);

  //  nc_close(AAdata.nc_id_MCD);
  fclose(AAdata.outputfile);  
  fclose(Zdata.trjfileZ);  
  fclose(Zdata.trjfilThetaAA); 

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename TACCMfilename parmfilename outputfilenameAA trjfileZ trjfileThetaZ\n",progname);
}
