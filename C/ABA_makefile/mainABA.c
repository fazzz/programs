#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "ABA.h"

#include "PTL.h"
#include "FFL.h"
#include "EF.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);
double dist(double *q,double *qgoal,int numclut);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d;
  int numatom,numclut,interval=1000,intervalflag;
  double dt=0.001;
  CLT *clt;
  double *Q,*frc,pot;
  double *q,*qgoal;
  double *qacc,*qvel,*qrot;
  double *predict,*correct;
  double kc=100.0,kd=50.0,ea=0.0;
  double KE,KEv,PEv;

  int MODE=NVT;
  double Tobj=10,KEobj;
  double omega=100.0,s_NVT=1.0,q_NVT=1.0,qvel_NVT=0.0,qacc_NVT=0.0,*predict_NVT,*correct_NVT;
  double k_B=1.98723e-3;

  double GoalReachedThershold=0.01,distance;

  int fflag=OFF;
  int dfalg=OFF,efalg=OFF,LJfalg=OFF,e14falg=OFF,LJ14falg=OFF;
  double coff=1.0;
  int *numclutparent,*terminal,*origin;

  double /***crd_nc,*/*energy;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;

  int *pairs;
  double *crd,*crdgoal,*mass;
  struct potential e;
  struct force f;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilenamestat,*inputfilenamegoal,*outputfilename,*parmfilename,*clustfilename;
  char *trjfilename;
  FILE *inputfilestat,*inputfilegoal,*outputfile,*parmfile,*clustfile;
  FILE *outputfile2;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"f",0,NULL,'f'},
    {"nve",0,NULL,'*'},
    {"dih",0,NULL,'d'},
    {"els",0,NULL,'e'},
    {"lj",0,NULL,'l'},
    {"e14",0,NULL,'1'},
    {"l14",0,NULL,'4'},
    {"h",0,NULL,'h'},
    {"kc",1,NULL,'k'},
    {"cof",1,NULL,'c'},
    {"temp",1,NULL,'t'},
    {"t",1,NULL,'t'},
    {"omg",1,NULL,'o'},
    {"int",1,NULL,'i'},
    {"GT",1,NULL,'g'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"f*del14ht:g:k:c:t:o:i:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case '*':
      MODE=NVE;
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
    case 'k':
      kc=atof(optarg);
      break;
    case 't':
      Tobj=atof(optarg);
      break;
    case 'o':
      omega=atof(optarg);
      break;
    case 'i':
      interval=atoi(optarg);
      break;
    case 'g':
      GoalReachedThershold=atof(optarg);
      break;
    case 'c':
      coff=atof(optarg);
      break;
    case 'f':
      fflag=ON;
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
  inputfilenamestat = *argv;
  inputfilenamegoal = *++argv;
  clustfilename     = *++argv;
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
  inputfilestat=efopen(inputfilenamestat,"r");
  getline(&line,&len,inputfilestat);
  fscanf(inputfilestat,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfilestat,"%lf",&crd[i*3+j]);
  fclose(inputfilestat);

  crdgoal=(double *)gcemalloc(sizeof(double)*numatom*3);
  inputfilegoal=efopen(inputfilenamegoal,"r");
  getline(&line,&len,inputfilegoal);
  fscanf(inputfilegoal,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfilegoal,"%lf",&crdgoal[i*3+j]);
  fclose(inputfilegoal);

  clustfile=efopen(clustfilename,"r");
  clt=ABAp_clustscan(clustfile,&numclut);
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

  ABAs_local_reference(clt,numclut,numatom,crd);
  ABAs_trans_Matrix(clt,numclut,numatom,crd);
  ABAs_inertia_matrix(clt,numclut,numatom,crd,mass);

  Q=(double *)gcemalloc(sizeof(double)*numclut);
  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  ffL_set_calcffandforce(&e,&f);

  q=(double *)gcemalloc(sizeof(double)*numclut);
  qgoal=(double *)gcemalloc(sizeof(double)*numclut);
  qacc=(double *)gcemalloc(sizeof(double)*numclut);
  qvel=(double *)gcemalloc(sizeof(double)*numclut);
  qrot=(double *)gcemalloc(sizeof(double)*numclut);
  predict=(double *)gcemalloc(sizeof(double)*numclut*6*6);
  correct=(double *)gcemalloc(sizeof(double)*numclut*6*6);
  predict_NVT=(double *)gcemalloc(sizeof(double)*numclut*6);
  correct_NVT=(double *)gcemalloc(sizeof(double)*numclut*6);

  for (i=0;i<numclut;++i) qvel[i]=0.00/*0.01*/;

  ABA_integ_set(q,qvel,predict,correct,numclut,dt);
  ABA_integ_set_NVT(q_NVT,qvel_NVT,predict_NVT,correct_NVT,dt);
  ABA_set_ini(q,clt,crd,numclut);
  ABA_set_ini(qgoal,clt,crdgoal,numclut);

  myncL_create_def_MCD(trjfilename,numatom,&nc_id_MCD);
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen("q.txt","w");
  for (i=0;(distance=dist(q,qgoal,numclut))>GoalReachedThershold && i<100000;++i) {
    ABA_integ_pret(qrot,qvel,q,predict,correct,dt,numclut);
    ABA_integ_pret_NVT(&qvel_NVT,&q_NVT,predict_NVT,correct_NVT,dt);
    ABA_update(clt,crd,qrot,numclut,numatom);

    for (j=0;j<numclut;++j) Q[j]=0.0;
    intervalflag=i%interval;
    ea=ABA_calcAttTorque_H(Q,clt,crd,q,qgoal,kc,numclut,outputfile2,intervalflag);
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

    ABA_calcKineE(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,q_NVT,qvel_NVT,s_NVT,numclut,numatom,MODE);
    solverABA(qacc,qvel,clt,Q,frc,crd,numclut,numatom,q,q_NVT,qvel_NVT,&qacc_NVT,s_NVT,KE,KEobj,MODE);

    ABA_integ_cort(qrot,qvel,q,qacc,predict,correct,dt,numclut);
    ABA_integ_cort_NVT(&qvel_NVT,&q_NVT,qacc_NVT,predict_NVT,correct_NVT,dt);
    ABA_update(clt,crd,qrot,numclut,numatom);

    ABA_calcKineE(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,q_NVT,qvel_NVT,s_NVT,numclut,numatom,MODE);
    if (i%interval==0) {
      fprintf(outputfile,"%lf %lf %lf %lf %lf %lf \n",pot,KE,KEv,PEv,pot+KE+PEv+KEv,ea);
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crd[j*3+k];
      myncL_put_crd_ene_MCD(nc_id_MCD,l,crd_nc,e,ea);
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
  printf("%s [-h] inputfilenamestat inputfilenamegoal clustfilename parmfilename outputfilename\n",progname);
}

double dist(double *q,double *qgoal,int numclut) {
  int i;
  double dist=0.0,d;
  double pi;

  pi=acos(-1.0);

  for (i=0;i<numclut;++i) {
    d=q[i]-qgoal[i];
    if (fabs(d)>pi) {
      if (d>0.0) d=-(2.0*pi-fabs(d));
      else d=2.0*pi-fabs(d);
    }
    dist+=d*d;
  }
  dist=sqrt(dist)/numclut;

  return dist;
}
