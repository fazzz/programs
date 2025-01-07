#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "ABA.h"
#include "GOLMAA_hybrid.h"
#include "GOLMAA_hybrid_set.h"

#include "PTL.h"
#include "EF.h"
#include "NC.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d;
  int numatom,numres,numclut,numstep=100000;
  int interval=1000,intervalout,intervalnc=1000,intervalflag;
  int numstepspe,numsflag=OFF;
  int addflag=OFF;
  int **nb_matrix;
  double dt=0.001;
  CLT *clt;
  double *Q,*frc,pot;
  double *q;
  double *qacc,*qvel,*qrot;
  double *predict,*correct;
  double KE,KEv,PEv;
  double dummy;

  int MODE=NVT,MODEV=OFF,TERMMODE=OFF;
  double Tobj=300,KEobj;
  double k_B=1.98723e-3;

  double s=1.0,s_vel=0.0,s_acc,gzi,gzi_vel,predict_s[6],correct_s[6];
  double Q_NH,tau=0.1,tau2;
  double T;
  int DOF;

  double constant=1.0;
  
  int num_NC,*ncmap,*ncmap_aa;
  int *indexncb;
  double Q_NC;

  double coff=1.0;
  int *numclutparent,*terminal,*origin;

  double *delta_Term,*vel_Term,*acc_Term,*acc_Term2,**predict_Term,**predict_Term2,**correct_Term,**correct_Term2;

  double /***crd_nc,*/*energy;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_AMBER nc_id;

  double *crd,*refcrd,*mass;
  struct potential_GOLMAA e;
  double pot_d;
  double R_C_D=1.0;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*reffilename,*outputfilename,*parmfilename,*clustfilename;
  char *trjfilename;
  char *rstfilename="rstfile",*rstvelfilename="rstvelfile";
  FILE *inputfile,*reffile,*outputfile,*parmfile,*clustfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"h",0,NULL,'h'},
    {"a",0,NULL,'a'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"ha",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 's':
      numstepspe=atof(optarg);
      numsflag=ON;
      break;
    case 'a':
      addflag=ON;
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
  inputfilename     = *argv;
  reffilename       = *++argv;
  parmfilename      = *++argv;
  outputfilename    = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  numres=AP.NRES;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);
  reffile=efopen(reffilename,"r");
  getline(&line,&len,reffile);
  fscanf(reffile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(reffile,"%lf",&refcrd[i*3+j]);
  fclose(reffile);

  numstep=mync_get_present_step_AMBER(inputfilename,&nc_id);

  GOLMAA_hybrid_ff_set_calcff(&e_GOLM,refcrd,numatom);

  if ( addflag==OFF )
    outputfile=efopen(outputfilename,"w");
  else
    outputfile=efopen(outputfilename,"a");

  for (i=0;i<numstep;++i) {
    mync_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id,crd_nc);
    for (j=0;j<numatom;++j)for (k=0;k<3;++k)crd[j*3+k]=crd_nc[j][k];

    GOLMAAff_calcff(crd,numatom,&e,ON);
    pot=e.p_t;

    fprintf(outputfile,"%10.8lf \n",pot/k_B/Tobj);

  }
  fclose(outputfile);
  nc_close((nc_id.ncid));

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfilename parmfilename outputfilename\n",progname);
}


