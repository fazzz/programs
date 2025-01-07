#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "REMDCGAA_TACCM_MPI_2_test.h"
#include "REMD_functions.h"
#include "EXT_SYS.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;

  struct AADataforREMD_test AAdata, CGdata;
  struct AACGCommonDataforREMD_test Cdata;
  struct TACCMDataforREMD_test Zdata;
  double parameterCG=0.2;

  int addflag=OFF,nobetaflag=OFF,AMBERMODEflag=OFF;

  double T=300,T_sim=300;
  double k_B=1.98723e-3;
  double beta=1.0;

  double p_t=0.0;

  double KZ=10.0;

  int *indexTACCM,**pairsZ;
  double *theta,*f2;
  char *TACCMfilename,*Ztrjfilename;
  FILE *TACCMfile,*Ztrjfile;

  double x[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double pi;

  char *inputfilename,*inputfilename2,*reffilename,*outputfilename,*parmfilename;
  FILE *Thetafile,*reffile,*outputfile,*parmfile;

  char *progname;

  int opt_idx=1;

  pi=acos(-1.0);

  struct option long_opt[] = {
    {"Amber",0,NULL,'A'},
    {"tempobj",1,NULL,'t'},
    {"tempsim",1,NULL,'s'},
    {"a",0,NULL,'a'},
    {"nobeta",0,NULL,'n'},
    {"KZ",1,NULL,'K'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hnaAt:K:s:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'A':
      AMBERMODEflag=ON;
      break;
    case 't':
      T=atof(optarg);
      break;
    case 's':
      T_sim=atof(optarg);
      break;
    case 'a':
      addflag=ON;
      break;
    case 'n':
      nobetaflag=ON;
      break;
    case 'K':
      KZ=atof(optarg);
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

  argc-=optind;
  argv+=optind;

  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }
  inputfilename     = *argv;
  inputfilename2    = *argv;
  Ztrjfilename      = *++argv;
  parmfilename      = *++argv;
  TACCMfilename     = *++argv;
  outputfilename    = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  numres=AP.NRES;

  extend_test_system(numEXT);
  for (i=0;i<AP.NATOM;++i) AP.CHRG[i]=0.0;
  for (i=0;i<(AP.NTYPES)*(AP.NTYPES+1)/2;++i){
    AP.CN1[i]=0.0; AP.CN2[i]=0.0;
  }

  Cdata.numatom=AP.NATOM;
  j=0;  for (i=0;i<Cdata.numatom;++i)  if (strncmp(AP.IGRAPH[i],"H",1)==0)  ++j;
  Cdata.numheavyatom=Cdata.numatom-j;  Cdata.numres=AP.NRES;

  if (AMBERMODEflag==ON) Cdata.numstep=mync_get_present_step_AMBER(inputfilename,&nc_id);
  else Cdata.numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MD);

  TACCMfile=efopen(TACCMfilename,"r");
  fscanf(TACCMfile,"%d",&(Zdata.numZ));
  Zdata.pairs=(int **)gcemalloc(sizeof(int *)*(Zdata.numZ));
  for (i=0;i<Zdata.numZ;++i) Zdata.pairs[i]=(int *)gcemalloc(sizeof(int)*5);
  for (i=0;i<Zdata.numZ;++i) {
    for (j=0;j<4;++j) 
      fscanf(TACCMfile,"%d",&(Zdata.pairs[i][j]));
    fscanf(TACCMfile,"%d",&(Zdata.pairs[i][j]));
  }
  fclose(TACCMfile);
  Zdata.Z=(double *)gcemalloc(sizeof(double)*numZ);

  AAdata.crd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  CGdata.crd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);

  ffL_set_calcffandforce(&e,&f);

  Thetafile=efopen(Thetafilename,"r");
  Ztrjfile=efopen(Ztrjfilename,"r");
  if ( nobetaflag==OFF )
    beta=1.0/k_B*(1.0/T_sim-1.0/T);
  if ( addflag==OFF ) outputfile=efopen(outputfilename,"w");
  else if ( addflag==ON )  outputfile=efopen(outputfilename,"a");
  for (i=0;i<numstep;++i) {
    for (j=0;j<numZ;++j) fscanf(Ztrjfile,"%lf",&(Zdata.Z[j]));

    CE_TACCM_CGAA(AAData.crd,CGData.crd,ZData.Z,
		  (CData.numatom),ZData.numZ,
		  KZAA,KZCG,ZData.pairs,
		  &EAAn_Xi,&ECGn_Xi,&EZn_Xi);

    fprintf(outputfile,"%d %e \n",i,EZn_Xi*beta);
  }
  fclose(outputfile);
  fclose(Ztrjfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] Thetafilename Ztrjfilename parmfilename TACCMfilename outputfilename\n",progname);
}


