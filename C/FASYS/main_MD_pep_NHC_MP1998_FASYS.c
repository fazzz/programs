
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

  double parameterCG=0.2;

  double avePE,aveKE,aveT,varPE,varKE,varT;

  int nc=1;                          
  double T0AA;
  double k_B=1.98723e-3;             
  double UNITT=418.4070;             
  double KTAA,tau=0.01,tau2,pi;                    

  double dt,dt2,wdt2[3],wdt4[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*parmfilename;
  char *outputfilename,*trjfilename;

  FILE *inputfile,*parmfile,*outputfile;

  char *progname;
  int opt_idx=1;

  pi=acos(-1.0);

  struct option long_opt[] = {
    {"vMode",1,NULL,'v'},
    {"nums",1,NULL,'s'},
    {"int",1,NULL,'i'},
    {"numEXT",1,NULL,'E'},
    {"tau",1,NULL,'a'}, 
    {"dt",1,NULL,'x'},
    {"TAA",1,NULL,'T'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hs:ia:E:x:T",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 's':
      Cdata.numstep=atoi(optarg); break;
    case 'i':
      Cdata.interval=atoi(optarg);  break;
    case 'T':
      T0AA=atof(optarg);  break;
    case 'v':
      vMode=ON;        break;
    case 'E':
      numEXT=atoi(optarg);  break;
    case 'a':
      tau=atof(optarg); break;
    case 'x':
      dt=atof(optarg);   break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);   exit(1);
    }
  }

  progname=*argv;  argc-=optind;  argv+=optind;

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  parmfilename   = *++argv;
  outputfilename = *++argv;
  trjfilename    = *++argv;

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

  Cdata.mass=(double *)gcemalloc(sizeof(double)*Cdata.numatom);
  for (i=0;i<Cdata.numatom;++i) Cdata.mass[i]=AP.AMASS[i];

  AAdata.crd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  AAdata.vel=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);  fscanf(inputfile,"%d",&d);
  for (i=0;i<Cdata.numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&(AAdata.crd[i*3+j]));
  fclose(inputfile);

  for (i=1;i<numEXT;++i) 
    for (j=0;j<Cdata.numatom/numEXT*3;++j) 
      AAdata.crd[i*Cdata.numatom/numEXT*3+j]=AAdata.crd[j];

  if ( vMode==OFF ) {
    MD_Generate_inivelo(AAdata.vel,Cdata.mass,Cdata.numatom,k_B*T0AA*UNITT);
    AAdata.zeta=0.0;      AAdata.V_zeta=0.0;
  }  

  tau=tau/2.0/pi;  tau2=tau*tau;  
  KTAA=k_B*T0AA;  AAdata.NfKT=(3.0*Cdata.numatom+1)*KTAA*UNITT;  AAdata.Q=tau2*KTAA*UNITT*(3.0*Cdata.numatom);

  ffL_set_calcffandforce(&(AAdata.e),&(AAdata.f));
  AAdata.e.parm.numnb=0;  AAdata.e.parm.num14=0;

  MD_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);

  myncL_create_def_MCD(trjfilename,Cdata.numatom,&(AAdata.nc_id_MCD));
  AAdata.outputfile=efopen(outputfilename,"w");

  runMD_NHC_MP1998_Amber_AAFF_test(AAdata.crd,AAdata.vel,Cdata.mass,(Cdata.numatom),
				   &(AAdata.zeta),&(AAdata.V_zeta),AAdata.Q,
				   AAdata.e,AAdata.f,
				   AAdata.T,AAdata.NfKT,
				   Cdata.numstep,Cdata.interval,&l,
				   dt,dt2,wdt2,wdt4,nc,
				   &avePE,&aveKE,&aveT,
				   &varPE,&varKE,&varT,
				   UNITT,k_B,
				   AAdata.nc_id_MCD,AAdata.outputfile);
  
  //  nc_close(AAdata.nc_id_MCD);
  fclose(AAdata.outputfile);  
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parmfilename outputfilename trjfilename\n",progname);
}
