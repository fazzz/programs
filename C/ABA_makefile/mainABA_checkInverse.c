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

#define ON 1
#define OFF 0

void usage(char *progname);
double dist(double *q,double *qgoal,int numclut);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d;
  int numatom,numclut;
  CLT *clt;
  double *Q,*frc,pot;
  double *q;
  double *qacc,*qvel,*qrot;
  double dummy;

  int dfalg=OFF,efalg=OFF,LJfalg=OFF,e14falg=OFF,LJ14falg=OFF;
  int *numclutparent,*terminal,*origin;

  int *pairs;
  double *crd,*crdgoal,*mass;
  struct potential e;
  struct force f;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilenamestat,*inputfilenamegoal,*parmfilename,*clustfilename;
  char *trjfilename;
  FILE *inputfilestat,*inputfilegoal,*parmfile,*clustfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"dih",0,NULL,'d'},
    {"els",0,NULL,'e'},
    {"lj",0,NULL,'l'},
    {"e14",0,NULL,'1'},
    {"l14",0,NULL,'4'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"del14h",long_opt,&opt_idx))!=-1) {
    switch(c) {
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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  inputfilenamestat = *argv;
  clustfilename     = *++argv;
  parmfilename      = *++argv;

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

  ABAs_local_reference(clt,numclut,numatom,crd);
  ABAs_trans_Matrix(clt,numclut,numatom,crd);
  ABAs_inertia_matrix(clt,numclut,numatom,crd,mass);

  Q=(double *)gcemalloc(sizeof(double)*numclut);
  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  ffL_set_calcffandforce(&e,&f);

  q=(double *)gcemalloc(sizeof(double)*numclut);
  qacc=(double *)gcemalloc(sizeof(double)*numclut);
  qvel=(double *)gcemalloc(sizeof(double)*numclut);

  for (i=0;i<numclut;++i) qvel[i]=0.01;

  for (j=0;j<numclut;++j) Q[j]=0.0;
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
      for (k=0;k<3;++k) frc[j*3+k]+=f.f_e[j*3+k];
    }
    if (LJfalg  == ON) {
      for (k=0;k<3;++k) frc[j*3+k]+=f.f_LJ[j*3+k];
    }
    if (e14falg == ON) {
      for (k=0;k<3;++k) frc[j*3+k]+=f.f_e_14[j*3+k];
    }
    if (LJ14falg== ON) {
      for (k=0;k<3;++k) frc[j*3+k]+=f.f_LJ_14[j*3+k];
    }
  }
  if (efalg   == ON)  pot+=0.5*e.p_e_t;
  if (LJfalg  == ON)  pot+=0.5*e.p_LJ_t;
  if (e14falg == ON)  pot+=0.5*e.p_e_14_t;
  if (LJ14falg== ON)  pot+=0.5*e.p_LJ_14_t;

  solverABA(qacc,qvel,clt,Q,frc,crd,numclut,numatom,q,0,0,&dummy,0,0,0,NVE);

  solverABA_Inverse(qacc,qvel,clt,Q,frc,crd,numclut,numatom);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfilename parmfilename \n",progname);
}


