#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PTL.h"
#include "FFL.h"
#include "EF.h"
#include "RAND.h"
#include "BOXMULL.h"
#include "MD.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

#define NVT 1
#define NVE 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d;
  int numatom,numstep=10000,interval=100;
  double dt=0.001;
  double *frc,PE;

  int bflag=ON,aflag=ON,dflag=ON,eflag=ON,LJflag=ON,e14flag=ON,LJ14flag=ON,natflag=ON,nnatflag=ON;

  int MODE=NVT;
  double Tobj=300,K,KE;
  double k_B=1.98723e-3,IsoCoff;
  double UNITT=418.4070;

  double *energy;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;

  double *crd,*mass,*vel,*vel0;

  struct potential e;
  struct force f;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*outputfilename,*parmfilename;
  char *trjfilename;
  FILE *inputfilestat,*outputfile,*parmfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"nve",0,NULL,'*'},
    {"bnd",0,NULL,'b'},
    {"ang",0,NULL,'a'},
    {"dih",0,NULL,'d'},
    {"e14",0,NULL,'1'},
    {"l14",0,NULL,'4'},
    {"e",0,NULL,'e'},
    {"l",0,NULL,'l'},
    {"h",0,NULL,'h'},
    {"nums",1,NULL,'s'},
    {"temp",1,NULL,'t'},
    {"t",1,NULL,'t'},
    {"int",1,NULL,'i'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"bad14el*+hs:t:i:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case '*':
      MODE=NVE;
      break;
    case 'b':
      bflag=OFF;
      break;
    case 'a':
      aflag=OFF;
      break;
    case 'd':
      dflag=OFF;
      break;
    case '1':
      e14flag=OFF;
      break;
    case '4':
      LJ14flag=OFF;
      break;
    case 'e':
      eflag=OFF;
      break;
    case 'l':
      LJflag=OFF;
      break;
    case 's':
      numstep=atoi(optarg);
      break;
    case 't':
      Tobj=atof(optarg);
      break;
    case 'i':
      interval=atoi(optarg);
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  parmfilename      = *++argv;
  outputfilename    = *++argv;
  trjfilename       = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  vel=(double *)gcemalloc(sizeof(double)*numatom*3);
  vel0=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfilestat=efopen(inputfilename,"r");
  getline(&line,&len,inputfilestat);
  fscanf(inputfilestat,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfilestat,"%lf",&crd[i*3+j]);
  fclose(inputfilestat);

  K=0.5*(3*numatom-1.0)*k_B*Tobj*UNITT;
  KE=MD_Generate_inivelo(vel,mass,numatom,k_B*Tobj*UNITT);
  for (i=0;i<numatom*3;++i) vel0[i]=vel[i];

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  ffL_set_calcffandforce(&e,&f);
  ffL_calcffandforce(crd,numatom,&e,&f);

  myncL_create_def_MCD(trjfilename,numatom,&nc_id_MCD);
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    KE=MD_Propagetor_vV_NVE_wc(crd,vel,mass,numatom,dt,&e,&f,bflag,aflag,dflag,eflag,LJflag,e14flag,LJ14flag);

    if (i%interval==0) {
      PE=0.0;
      if(bflag==ON)  PE+=e.p_b_t;
      if(aflag==ON)  PE+=e.p_a_t;
      if(dflag==ON)  PE+=e.p_d_t;
      if (eflag == ON)   PE+=0.5*e.p_e_t;
      if (LJflag== ON)   PE+=0.5*e.p_LJ_t;
      if (e14flag == ON) PE+=0.5*e.p_e_14_t;
      if (LJ14flag== ON) PE+=0.5*e.p_LJ_14_t;
	
      fprintf(outputfile,"%d %e %e %e %e\n",i+1,PE,KE/UNITT,PE+KE/UNITT,KE/UNITT/(0.5*numatom*3*k_B));
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crd[j*3+k];
      myncL_put_crd_ene_MCD(nc_id_MCD,l,crd_nc,e,0.0);
      ++l;
    }
  }
  fclose(outputfile);
  nc_close((nc_id_MCD.ncid));

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parmfilename outputfilename\n",progname);
}


