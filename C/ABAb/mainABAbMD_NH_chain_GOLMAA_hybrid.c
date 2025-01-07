#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "ABAb.h"
#include "GOLMAA_hybrid_set.h"
#include "GOLMAA_hybrid.h"

#include "PTL.h"
#include "FFL.h"
#include "EF.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d;
  int numatom,numres,numclut,numstep=100000;
  int interval=1000,intervalout=1000,intervalnc=1000,intervalflag;
  double dt=0.001;
  CLTb *clt;
  double *Q,*frc,pot;
  double *q;
  double *qacc,*qvel,*qrot;
  double *predict,*correct;
  double KE,KEv,PEv;
  double dummy;

  int MODE=NVT,MODEV=OFF,TERMMODE=OFF;
  double Tobj=300,KEobj,KBT;
  double k_B=1.98723e-3;

  double *zeta,*zeta_vel,*zeta_acc,**predict_zeta,**correct_zeta;
  double *Q_NH,tau=0.1,tau2;
  double T;
  int DOF,M=4;

  int dfalg=ON,efalg=ON,LJfalg=ON,e14falg=ON,LJ14falg=ON,natfalg=ON,nnatfalg=ON;
  double coff=1.0;
  int *numclutparent,*terminal,*origin;

  double *delta_Term,*vel_Term,*acc_Term,*acc_Term2,**predict_Term,**predict_Term2,**correct_Term,**correct_Term2;

  double /***crd_nc,*/*energy;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_AMBER nc_id;

  double *crd,*refcrd,*mass;
  struct potential e;
  struct force f;
  struct potential_GOLMAA_hybrid e_GOLM;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*refcrdfilename;
  char *outputfilename,*outputfilename2,*parmfilename,*clustfilename,*inputvelofilename;
  char *trjfilename;
  char *rstfilename="rstfile",*rstvelfilename="rstvelfile";
  FILE *inputfile,*reffile,*inputfilegoal,*outputfile,*outputfile2,*parmfile,*clustfile,*inputvelofile;
  FILE *velfile,*rstfile,*rstvelfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"f",0,NULL,'f'},
    {"nve",0,NULL,'*'},
    {"termon",0,NULL,'+'},
    {"dih",0,NULL,'d'},
    {"e14",0,NULL,'1'},
    {"l14",0,NULL,'4'},
    {"nat",0,NULL,'c'},
    {"repul",0,NULL,'N'},
    {"h",0,NULL,'h'},
    {"kc",1,NULL,'k'},
    {"cof",1,NULL,'c'},
    {"temp",1,NULL,'t'},
    {"dt",1,NULL,'@'},
    {"omg",1,NULL,'o'},
    {"nums",1,NULL,'s'},
    {"int",1,NULL,'i'},
    {"intnc",1,NULL,'j'},
    {"intout",1,NULL,'k'},
    {"vel",1,NULL,'v'},
    {"tau",1,NULL,'?'},
    {"nchain",1,NULL,'M'},
    {"rst",1,NULL,'{'},
    {"rstv",1,NULL,'}'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"*+d14cNh@:t:k:c:t:o:s:i:j:k:v:?:M:{:}:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case '*':
      MODE=NVE;
      break;
    case '+':
      TERMMODE=ON;
      break;
    case 'd':
      dfalg=OFF;
      break;
    case '1':
      e14falg=OFF;
      break;
    case '4':
      LJ14falg=OFF;
      break;
    case 'c':
      natfalg=OFF;
      break;
    case 'N':
      nnatfalg=OFF;
      break;
    case 't':
      Tobj=atof(optarg);
      break;
    case 'o':
      tau=atof(optarg);
      break;
    case '@':
      dt=atof(optarg);
      break;
    case 'i':
      interval=atoi(optarg);
      break;
    case 'j':
      intervalnc=atoi(optarg);
      break;
    case 'k':
      intervalout=atoi(optarg);
      break;
    case 's':
      numstep=atoi(optarg);
      break;
    case 'v':
      inputvelofilename=optarg;
      MODEV=ON;
      break;
    case '?':
      tau=atof(optarg);
      break;
    case 'M':
      M=atoi(optarg);
      break;
    case '{':
      rstfilename=optarg;
      break;
    case '}':
      rstvelfilename=optarg;
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

  if (argc < 7) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  refcrdfilename = *++argv;
  clustfilename  = *++argv;
  parmfilename   = *++argv;
  outputfilename = *++argv;
  outputfilename2= *++argv;
  trjfilename    = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  numres=AP.NRES;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);
  reffile=efopen(refcrdfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(reffile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(reffile,"%lf",&refcrd[i*3+j]);
  fclose(reffile);

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

  ABAbs_local_reference(clt,numclut,numatom,crd);
  ABAbs_trans_Matrix(clt,numclut,numatom,crd);
  ABAbs_inertia_matrix(clt,numclut,numatom,crd,mass);

  Q=(double *)gcemalloc(sizeof(double)*numclut);
  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  ffL_set_calcffandforce(&e,&f);

  //  GOLMAA_hybrid_ff_set_calcff(&e_GOLM,refcrd,numatom,numres,e.parm.indexnb,e.parm.numnb);
  GOLMAA_hybrid_ff_set_calcff_3(&e_GOLM,refcrd,numatom,numres,e.parm.indexnb,e.parm.numnb);

  q=(double *)gcemalloc(sizeof(double)*numclut);
  qacc=(double *)gcemalloc(sizeof(double)*numclut);
  qvel=(double *)gcemalloc(sizeof(double)*numclut);
  qrot=(double *)gcemalloc(sizeof(double)*numclut);
  predict=(double *)gcemalloc(sizeof(double)*numclut*6/**6*/);
  correct=(double *)gcemalloc(sizeof(double)*numclut*6/**6*/);

  zeta=(double *)gcemalloc(sizeof(double)*M);
  zeta_vel=(double *)gcemalloc(sizeof(double)*M);
  zeta_acc=(double *)gcemalloc(sizeof(double)*M);
  predict_zeta=(double **)gcemalloc(sizeof(double *)*M);
  correct_zeta=(double **)gcemalloc(sizeof(double *)*M);
  for (i=0;i<M;++i) {
    predict_zeta[i]=(double *)gcemalloc(sizeof(double)*6);
    correct_zeta[i]=(double *)gcemalloc(sizeof(double)*6);
  }
  Q_NH=(double *)gcemalloc(sizeof(double)*M);

  delta_Term=(double *)gcemalloc(sizeof(double)*6);
  vel_Term=(double *)gcemalloc(sizeof(double)*6);
  acc_Term=(double *)gcemalloc(sizeof(double)*6);
  acc_Term2=(double *)gcemalloc(sizeof(double)*6);
  predict_Term=(double **)gcemalloc(sizeof(double *)*6);
  predict_Term2=(double **)gcemalloc(sizeof(double *)*6);
  correct_Term=(double **)gcemalloc(sizeof(double *)*6);
  correct_Term2=(double **)gcemalloc(sizeof(double *)*6);
  for (i=0;i<6;++i) {
    predict_Term[i]=(double *)gcemalloc(sizeof(double)*6);
    predict_Term2[i]=(double *)gcemalloc(sizeof(double)*6);
    correct_Term[i]=(double *)gcemalloc(sizeof(double)*6);
    correct_Term2[i]=(double *)gcemalloc(sizeof(double)*6);
  }

  for (i=0;i<6;++i) vel_Term[i]=0.0;
  for (i=0;i<numclut;++i) qvel[i]=0.000/*1.00*/;

  DOF=(numclut-1);
  if (TERMMODE==ON) DOF+=6;
  KEobj=0.5*DOF*k_B*Tobj;
  KBT=k_B*Tobj;

  ABAbNH_chain_set_new(tau,&tau2,Q_NH,M,DOF,KBT,correct_zeta);

  ABAb_integ_set(q,qvel,predict,correct,numclut,dt);

  if (MODEV==ON) 
    ABAbs_restat_read_chain(inputvelofilename,numclut,correct,correct_Term,correct_Term2,correct_zeta,M,MODE,TERMMODE);
  
  myncL_create_def_AMBER(trjfilename,numatom,&nc_id);
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");

  for (i=0;i<numstep;++i) {
    ABAb_integ_pret(qrot,qvel,q,predict,correct,dt,numclut);
    if (TERMMODE==ON) 
      ABAb_integ_pret_Term(predict_Term,predict_Term2,correct_Term,correct_Term2,vel_Term,delta_Term,dt);
    if (MODE==NVT) 
      ABAbNH_chain_update_pret(zeta,zeta_vel,predict_zeta,correct_zeta,M,dt);

    ABAb_update(clt,crd,qrot,numclut,numatom);
    //    if (TERMMODE==ON) ABAb_update_Term(crd,delta_Term,numatom);

    for (j=0;j<numclut;++j) Q[j]=0.0;
    for (j=0;j<numatom;++j) for (k=0;k<3;++k) frc[j*3+k]=0.0;
    intervalflag=i%interval;

    if(dfalg==ON) ffL_calcTorque_woH(Q,crd,numclut,numclutparent,terminal,origin);
    if (e14falg == ON || LJ14falg == ON ) {
      //      ffL_calcffandforce_14D_woH(crd,numatom,&e,&f);
      ffL_calcffandforce_14vdWDAB_woH(crd,numatom,&e,&f);
    }
    for (j=0;j<numatom;++j) {
      /**********************************************************************/
      /* if (e14falg == ON) {						    */
      /* 	for (k=0;k<3;++k) frc[j*3+k]+=/\*coff**\/f.f_e_14[j*3+k];   */
      /* }								    */
      /**********************************************************************/
      if (LJ14falg== ON) {
    	for (k=0;k<3;++k) frc[j*3+k]+=/*coff**/f.f_LJ_14[j*3+k];
      }
    }

    if (natfalg == ON || nnatfalg == ON ) 
      GOLMAA_hyb_ff_calcff(crd,numatom,&e_GOLM);
    for (j=0;j<numatom;++j) {
       if (natfalg == ON)
	 for (k=0;k<3;++k) frc[j*3+k]+=e_GOLM.f_natatt[j][k];
       if (nnatfalg == ON)
	 for (k=0;k<3;++k) frc[j*3+k]+=e_GOLM.f_repul[j][k];
    }

    if (TERMMODE==ON) 
      ABAb_calcKineE_TermOn_new(&KE,clt,crd,qvel,vel_Term,numclut,numatom);
    else 
      ABAb_calcKineE_new(&KE,clt,crd,qvel,numclut,numatom);
    T=KE/(DOF*k_B)*2.0;

    if (TERMMODE==ON) 
      solverABAb_TermOn_NH_chain(qacc,qvel,clt,Q,frc,crd,numclut,numatom,acc_Term,acc_Term2,vel_Term,zeta_vel,zeta_acc,M,DOF,KBT,Q_NH,tau2,T,Tobj);
    else 
      solverABAb_NH_chain(qacc,qvel,clt,Q,frc,crd,numclut,numatom,zeta_vel,zeta_acc,M,DOF,KBT,Q_NH,tau2,T,Tobj);

    ABAb_integ_cort(qrot,qvel,q,qacc,predict,correct,dt,numclut);
    if (TERMMODE==ON)
      ABAb_integ_cort_Term(predict_Term,predict_Term2,correct_Term,correct_Term2,acc_Term,acc_Term2,vel_Term,delta_Term,dt);
    if (MODE==NVT) ABAbNH_chain_update_cort(zeta,zeta_vel,zeta_acc,predict_zeta,correct_zeta,M,dt);

    ABAb_update(clt,crd,qrot,numclut,numatom);
    //    if (TERMMODE==ON) ABAb_update_Term(crd,delta_Term,numatom);

    if (TERMMODE==ON)
      ABAb_calcKineE_TermOn_new(&KE,clt,crd,qvel,vel_Term,numclut,numatom);
    else
      ABAb_calcKineE(&KE,clt,crd,qvel,numclut,numatom);
    if (MODE==NVT)
      ABAbNH_chain_calcKE_new(zeta,zeta_vel,Q_NH,M,DOF,KBT,&PEv,&KEv);

    if (i%interval==0) {
      pot=0.0;
      if (dfalg==ON) pot=e.p_d_t;
      //      if (e14falg == ON)  pot+=0.5*e.p_e_14_t;
      if (LJ14falg== ON)  pot+=0.5*e.p_LJ_14_t;
      if (natfalg == ON)  pot+=0.5*e_GOLM.p_natatt_t;
      if (nnatfalg == ON) pot+=0.5*e_GOLM.p_repul_t;
      fprintf(outputfile,"%d %e %e %e %e %e %e %e\n",i,pot,KE,KEv,PEv,pot+KE+PEv+KEv,zeta[0],T);
    }

    if (i%intervalout==0) {
      pot=0.0;
      if (dfalg==ON) pot=e.p_d_t;
      //      if (e14falg == ON)  pot+=0.5*e.p_e_14_t;
      if (LJ14falg== ON)  pot+=0.5*e.p_LJ_14_t;
      if (natfalg == ON)  pot+=e_GOLM.p_natatt_t;
      if (nnatfalg == ON) pot+=e_GOLM.p_repul_t;
      fprintf(outputfile2,"E_t    = %e \n",pot+KE+PEv+KEv);
      fprintf(outputfile2,"KE     = %e \n",KE);
      fprintf(outputfile2,"KEv    = %e \n",KEv);
      fprintf(outputfile2,"PEv    = %e \n",PEv);
      fprintf(outputfile2,"p_t    = %e \n",pot);
      fprintf(outputfile2,"p_NC   = %e \n",e_GOLM.p_natatt_t);
      fprintf(outputfile2,"p_NNC  = %e \n",e_GOLM.p_repul_t);
      //      fprintf(outputfile2,"p_14es = %e \n",e.p_e_14_t);
      fprintf(outputfile2,"p_14LJ = %e \n",e.p_LJ_14_t);
      fprintf(outputfile2,"p_dihe = %e \n",e.p_d_t);
    }

    if (i%intervalnc==0) {
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crd[j*3+k];
      myncL_put_crd_AMBER(nc_id,l,crd_nc);
      ++l;
    }
  }
  fclose(outputfile);
  fclose(outputfile2);
  nc_close((nc_id.ncid));

  rstfile=efopen(rstfilename,"w");
  fprintf(rstfile,"ACE\n ");
  fprintf(rstfile,"%d \n",numatom);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      fprintf(rstfile,"%10.8lf ",crd[i*3+j]);
    }
    fprintf(rstfile,"\n");
  }
  fclose(rstfile);

  ABAbs_restat_write_vel_chain(rstvelfilename,numclut,correct,correct_Term,correct_Term2,correct_zeta,M,MODE,TERMMODE);
  fclose(outputfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename refcrdfilename clustfilename parmfilename outputfilename\n",progname);
}


