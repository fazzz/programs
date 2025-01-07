#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "ABAb.h"

#include "PTL.h"
#include "FFL.h"
#include "EF.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d;
  int numatom,numclut,numstep=100000,interval=1000,intervalflag;
  double dt=0.001;
  CLTb *clt;
  double *Q,*frc,pot;
  double *q,*qgoal;
  double *qacc,*qvel,*qrot;
  double *qvel_h,*qvel_oh,*qvel_hu;
  double KE,KEv,PEv;
  double dummy;

  int MODE=NVE,MODEV=OFF,TERMMDOE=ON;
  double Tobj=10,KEobj;
  double omega=100.0,s_NVT=1.0,q_NVT=1.0,qvel_NVT=0.0,qacc_NVT=0.0,*predict_NVT,*correct_NVT;
  double k_B=1.98723e-3;

  int dfalg=ON,efalg=ON,LJfalg=ON,e14falg=ON,LJ14falg=ON;
  double coff=1.0;
  int *numclutparent,*terminal,*origin;

  double *vel_Term,*acc_Term,*acc_Term2,*vel_Term_h,*vel_Term_oh,*vel_Term_hu;

  double /***crd_nc,*/*energy;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;

  double *crd,*mass;
  struct potential e;
  struct force f;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilenamestat,*inputfilenamegoal,*outputfilename,*outputfilename2,*parmfilename,*clustfilename,*inputvelofilename;
  char *trjfilename;
  FILE *inputfilestat,*inputfilegoal,*outputfile,*outputfile2,*parmfile,*clustfile,*inputvelofile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"f",0,NULL,'f'},
    {"nve",0,NULL,'*'},
    {"termoff",0,NULL,'+'},
    {"dih",0,NULL,'d'},
    {"els",0,NULL,'e'},
    {"lj",0,NULL,'l'},
    {"e14",0,NULL,'1'},
    {"l14",0,NULL,'4'},
    {"h",0,NULL,'h'},
    {"kc",1,NULL,'k'},
    {"cof",1,NULL,'c'},
    {"temp",1,NULL,'t'},
    {"dt",1,NULL,'@'},
    {"omg",1,NULL,'o'},
    {"nums",1,NULL,'s'},
    {"int",1,NULL,'i'},
    {"vel",1,NULL,'v'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"*+del14h@:t:k:c:t:o:s:i:v:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case '*':
      MODE=NVE;
      break;
    case '+':
      TERMMDOE=OFF;
      break;
    case 'd':
      dfalg=ON;
      break;
    case 'e':
      efalg=ON;
      break;
    case 'l':
      LJfalg=ON;
      break;
    case '1':
      e14falg=ON;
      break;
    case '4':
      LJ14falg=ON;
      break;
    case 't':
      Tobj=atof(optarg);
      break;
    case 'o':
      omega=atof(optarg);
      break;
    case '@':
      dt=atof(optarg);
      break;
    case 'i':
      interval=atoi(optarg);
      break;
    case 's':
      numstep=atoi(optarg);
      break;
    case 'c':
      coff=atof(optarg);
      break;
    case 'v':
      inputvelofilename=optarg;
      MODEV=ON;
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

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  inputfilenamestat = *argv;
  clustfilename     = *++argv;
  parmfilename      = *++argv;
  outputfilename    = *++argv;
  outputfilename2    = *++argv;
  trjfilename       = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  inputfilestat=efopen(inputfilenamestat,"r");
  getline(&line,&len,inputfilestat);
  fscanf(inputfilestat,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfilestat,"%lf",&crd[i*3+j]);
  fclose(inputfilestat);

  clustfile=efopen(clustfilename,"r");
  clt=ABAbp_clustscan(clustfile,&numclut);
  fclose(clustfile);

  numclutparent=(int *)gcemalloc(sizeof(int)*numclut);
  terminal=(int *)gcemalloc(sizeof(int)*numclut);
  origin=(int *)gcemalloc(sizeof(int)*numclut);
  for (i=0;i<numclut;++i) {
    numclutparent[i]=clt[i].nNumClutOfParent;
    terminal[i]=clt[i].terminal_atom_a[0];
    origin[i]=clt[i].origin_atom_a;
  }

  KEobj=(numclut-1)*k_B*Tobj;
  s_NVT=(numclut-1)*k_B*Tobj/UNIT*omega*omega;

  ABAbs_local_reference(clt,numclut,numatom,crd);
  ABAbs_trans_Matrix(clt,numclut,numatom,crd);
  ABAbs_inertia_matrix(clt,numclut,numatom,crd,mass);

  Q=(double *)gcemalloc(sizeof(double)*numclut);
  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  ffL_set_calcffandforce(&e,&f);

  q=(double *)gcemalloc(sizeof(double)*numclut);
  qgoal=(double *)gcemalloc(sizeof(double)*numclut);
  qacc=(double *)gcemalloc(sizeof(double)*numclut);
  qvel=(double *)gcemalloc(sizeof(double)*numclut);
  qrot=(double *)gcemalloc(sizeof(double)*numclut);
  qvel_h=(double *)gcemalloc(sizeof(double)*numclut);
  qvel_oh=(double *)gcemalloc(sizeof(double)*numclut);
  qvel_hu=(double *)gcemalloc(sizeof(double)*numclut);

  vel_Term=(double *)gcemalloc(sizeof(double)*6);
  acc_Term=(double *)gcemalloc(sizeof(double)*6);
  acc_Term2=(double *)gcemalloc(sizeof(double)*6);
  vel_Term_h=(double *)gcemalloc(sizeof(double)*6);
  vel_Term_hu=(double *)gcemalloc(sizeof(double)*6);
  vel_Term_oh=(double *)gcemalloc(sizeof(double)*6);

  for (i=0;i<6;++i) {
    vel_Term[i]=0.0;
    vel_Term_h[i]=0.0;
    vel_Term_hu[i]=0.0;
    vel_Term_oh[i]=0.0;
  }
  for (i=0;i<numclut;++i) {
    qvel[i]=0.00/*0.01*/;
    qvel_h[i]=0.00/*0.01*/;
    qvel_hu[i]=0.00/*0.01*/;
    qvel_oh[i]=0.00/*0.01*/;
  }
  if (MODEV==ON) {
    inputvelofile=efopen(inputvelofilename,"r");
    for (i=0;i<6;++i) fscanf(inputvelofile,"%lf",&vel_Term[i]);
    for (i=1;i<numclut;++i) fscanf(inputvelofile,"%lf",&qvel[i]);
    fclose(inputvelofile);
  }

  //  myncL_create_def_MCD(trjfilename,numatom,&nc_id_MCD);
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<numstep;++i) {
    for (j=0;j<numclut;++j) qvel[j]=1.5*qvel_h[j]-0.5*qvel_oh[j];
    for (j=0;j<6;++j) vel_Term[j]=1.5*vel_Term_h[j]-0.5*vel_Term_oh[j];

    for (j=0;j<numclut;++j) Q[j]=0.0;
    intervalflag=i%interval;
    if (dfalg==ON) {
      ffL_calcTorque(Q,crd,numclut,numclutparent,terminal,origin);
    }
    ffL_calcffandforce(crd,numatom,&e,&f);
    pot=0.0;
    if (dfalg==ON) 
      pot=e.p_d_t;
    for (j=0;j<numatom;++j) for (k=0;k<3;++k) frc[j*3+k]=0.0;
    for (j=0;j<numatom;++j) {
      if (efalg   == ON) {
	for (k=0;k<3;++k) frc[j*3+k]+=coff*f.f_e[j*3+k];
      }
      if (LJfalg  == ON) {
	for (k=0;k<3;++k) frc[j*3+k]+=coff*f.f_LJ[j*3+k];
      }
      if (e14falg == ON) {
	for (k=0;k<3;++k) frc[j*3+k]+=coff*f.f_e_14[j*3+k];
      }
      if (LJ14falg== ON) {
	for (k=0;k<3;++k) frc[j*3+k]+=coff*f.f_LJ_14[j*3+k];
      }
    }
    if (efalg   == ON)  pot+=0.5*e.p_e_t;
    if (LJfalg  == ON)  pot+=0.5*e.p_LJ_t;
    if (e14falg == ON)  pot+=0.5*e.p_e_14_t;
    if (LJ14falg== ON)  pot+=0.5*e.p_LJ_14_t;

    ABAb_calcKineE(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,q_NVT,qvel_NVT,s_NVT,numclut,numatom,MODE);
    if (TERMMDOE==OFF)
      solverABAb(qacc,qvel,clt,Q,frc,crd,numclut,numatom,q,q_NVT,qvel_NVT,&qacc_NVT,s_NVT,KE,KEobj,MODE);
    if (TERMMDOE==ON)
      solverABAb_TermOn(qacc,qvel,clt,Q,frc,crd,numclut,numatom,q,q_NVT,qvel_NVT,&qacc_NVT,s_NVT,acc_Term,acc_Term2,vel_Term,KE,KEobj,MODE);

    for (j=0;j<numclut;++j) qvel_hu[j]=qvel_h[j]+dt*qacc[j];
    for (j=0;j<6;++j) vel_Term_hu[j]=vel_Term_h[j]+dt*acc_Term[j];

    for (j=0;j<numclut;++j) qvel[j]=0.5*qvel_hu[j]+0.5*qvel_h[j];
    for (j=0;j<6;++j) vel_Term[j]=0.5*vel_Term_hu[j]+0.5*vel_Term_h[j];

    for (j=0;j<numclut;++j) qrot[j]=dt*qvel_hu[j];

    ABAb_update(clt,crd,qrot,numclut,numatom);

    if (TERMMDOE==OFF)
      ABAb_calcKineE(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,q_NVT,qvel_NVT,s_NVT,numclut,numatom,MODE,numclut);
    if (TERMMDOE==ON)
      ABAb_calcKineE_TermOn(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,q_NVT,qvel_NVT,s_NVT,vel_Term,numclut,numatom,MODE,numclut+6);

    for (j=0;j<numclut;++j) {
      qvel_oh[j]=qvel_h[j];
      qvel_h[j]=qvel_hu[j];
    }
    for (j=0;j<6;++j) {
      vel_Term_oh[j]=vel_Term_h[j];
      vel_Term_h[j]=vel_Term_hu[j];
    }

    if (i%interval==0) {
      fprintf(outputfile,"%lf %lf %lf %lf %lf \n",pot,KE,KEv,PEv,pot+KE+PEv+KEv);
      ABAb_out_formated(outputfile2,pot,KE,i,dt);
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crd[j*3+k];
      //      myncL_put_crd_ene_MCD(nc_id_MCD,l,crd_nc,e,ea);
      ++l;
    }
  }
  fclose(outputfile);
  fclose(outputfile2);
  nc_close((nc_id_MCD.ncid));

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilenamestat clustfilename parmfilename outputfilename\n",progname);
}


