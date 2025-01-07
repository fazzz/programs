#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "ABA.h"
#include "GOLMAA.h"
#include "GOLMAA_set.h"

#include "PTL.h"
#include "EF.h"
#include "NC.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d;
  int numatom,numres,numclut,numstep=100000,interval=1000,intervalnc=1000,intervalflag;
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

  double *delta_Term,*vel_Term,*acc_Term,*acc_Term2,**predict_Term,**predict_Term2,**correctt_Term,**correct_Term2;

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

  char *inputfilename,*reffilename,*outputfilename,*outputfilename2,*parmfilename,*clustfilename,*inputvelofilename;
  char *trjfilename;
  char *rstfilename="rstfile",*rstvelfilename="rstvelfile";
  FILE *inputfile,*reffile,*inputfilegoal,*outputfile,*outputfile2,*parmfile,*clustfile,*inputvelofile;
  FILE *velfile,*rstfile,*rstvelfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"nve",0,NULL,'*'},
    {"termon",0,NULL,'+'},
    {"h",0,NULL,'h'},
    {"temp",1,NULL,'t'},
    {"dt",1,NULL,'@'},
    {"nums",1,NULL,'s'},
    {"int",1,NULL,'i'},
    {"intnc",1,NULL,'j'},
    {"vel",1,NULL,'v'},
    {"rate",1,NULL,'|'},
    {"tau",1,NULL,'?'},
    {"rst",1,NULL,'{'},
    {"rstv",1,NULL,'}'},
    {"cons",1,NULL,'>'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"*+h@:t:k:c:t:i:j:v:|:?:{:}:>:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case '*':
      MODE=NVE;
      break;
    case '+':
      TERMMODE=ON;
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
    case '|':
      R_C_D=atof(optarg);
      break;
    case '{':
      rstfilename=optarg;
      break;
    case '}':
      rstvelfilename=optarg;
      break;
    case '>':
      constant=atof(optarg);
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
  inputfilename     = *argv;
  reffilename       = *++argv;
  clustfilename     = *++argv;
  parmfilename      = *++argv;
  outputfilename    = *++argv;
  outputfilename2   = *++argv;
  trjfilename       = *++argv;

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
  reffile=efopen(reffilename,"r");
  getline(&line,&len,reffile);
  fscanf(reffile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(reffile,"%lf",&refcrd[i*3+j]);
  fclose(reffile);

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
  nb_matrix=(int **)gcemalloc(sizeof(int *)*numatom);
  GOLMAAff_set_calcff(&e,refcrd,numatom,nb_matrix,R_C_D,constant);

  q=(double *)gcemalloc(sizeof(double)*numclut);
  qacc=(double *)gcemalloc(sizeof(double)*numclut);
  qvel=(double *)gcemalloc(sizeof(double)*numclut);
  qrot=(double *)gcemalloc(sizeof(double)*numclut);
  predict=(double *)gcemalloc(sizeof(double)*numclut*6*6);
  correct=(double *)gcemalloc(sizeof(double)*numclut*6*6);

  delta_Term=(double *)gcemalloc(sizeof(double)*6);
  vel_Term=(double *)gcemalloc(sizeof(double)*6);
  acc_Term=(double *)gcemalloc(sizeof(double)*6);
  acc_Term2=(double *)gcemalloc(sizeof(double)*6);
  predict_Term=(double **)gcemalloc(sizeof(double *)*6);
  predict_Term2=(double **)gcemalloc(sizeof(double *)*6);
  correctt_Term=(double **)gcemalloc(sizeof(double *)*6);
  correct_Term2=(double **)gcemalloc(sizeof(double *)*6);
  for (i=0;i<6;++i) {
    predict_Term[i]=(double *)gcemalloc(sizeof(double)*6);
    predict_Term2[i]=(double *)gcemalloc(sizeof(double)*6);
    correctt_Term[i]=(double *)gcemalloc(sizeof(double)*6);
    correct_Term2[i]=(double *)gcemalloc(sizeof(double)*6);
  }

  for (i=0;i<6;++i) vel_Term[i]=0.0;
  for (i=0;i<numclut;++i) qvel[i]=0.00/*1.00*/;

  if (MODEV==ON) {
    inputvelofile=efopen(inputvelofilename,"r");
    if (TERMMODE==ON) for (i=0;i<6;++i) fscanf(inputvelofile,"%lf",&vel_Term[i]);
    for (i=0;i<numclut;++i) fscanf(inputvelofile,"%lf",&qvel[i]);
    if (MODE==NVT) fscanf(inputvelofile,"%lf %lf",&s,&s_vel);
    fclose(inputvelofile);
  }

  DOF=(numclut-1);
  if (TERMMODE==ON) DOF+=6;
  KEobj=/*0.5**/DOF*k_B*Tobj;
  ABANH_set(s,s_vel,gzi,predict_s,correct_s,tau,&tau2,&Q_NH,KEobj,dt);

  ABA_integ_set(q,qvel,predict,correct,numclut,dt);

  myncL_create_def_AMBER(trjfilename,numatom,&nc_id);
  outputfile=efopen(outputfilename,"w");
  //  outputfile2=efopen(outputfilename2,"w");
  ncmap_aa=(int *)gcemalloc(sizeof(int)*numatom*numatom);
  ncmap=(int *)gcemalloc(sizeof(int)*numres*numres);
  indexncb=make_native_contact_list(&num_NC,refcrd,numatom,numres,criteria_NC,ncmap);

  for (i=0;i<numstep;++i) {
    ABA_integ_pret(qrot,qvel,q,predict,correct,dt,numclut);
    if (TERMMODE==ON) ABA_integ_pret_Term(predict_Term,predict_Term2,correctt_Term,correct_Term2,vel_Term,delta_Term,dt);

    if (MODE==NVT)
      ABANH_update_pret(&s,&s_vel,predict_s,correct_s,dt);
    ABA_update(clt,crd,qrot,numclut,numatom);
    if (TERMMODE==ON) ABA_update_Term(crd,delta_Term,numatom);

    for (j=0;j<numclut;++j) Q[j]=0.0;
    intervalflag=i%interval;
    //    pot_d=GOLMAA_calcTorque(Q,crd,e.DEQ,e.FC_dihed,numclut,numclutparent,terminal,origin);
    GOLMAAff_calcff(crd,numatom,&e,OFF,ON,ON,nb_matrix);
    pot=e.p_t/*+pot_d*/;
    e.p_d_t=0.0/*pot_d*/;
    for (j=0;j<numatom;++j) for (k=0;k<3;++k) frc[j*3+k]=e.f_t[j][k];

    ABA_calcKineE_TermOn(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,s,s_vel,Q_NH,vel_Term,numclut,numatom,MODE,numclut+6);
    T=KE/(DOF*k_B)/**2.0*/;
    if (TERMMODE==ON) {
      solverABA_TermOn_NH(qacc,qvel,clt,Q,frc,crd,numclut,numatom,q,s,s_vel,&s_acc,tau2,acc_Term,acc_Term2,vel_Term,T,Tobj);
    }    
    else 
      solverABA_NH(qacc,qvel,clt,Q,frc,crd,numclut,numatom,q,s,s_vel,&s_acc,tau2,T,Tobj);
    ABA_integ_cort(qrot,qvel,q,qacc,predict,correct,dt,numclut);
    if (TERMMODE==ON)
      ABA_integ_cort_Term(predict_Term,predict_Term2,correctt_Term,correct_Term2,acc_Term,acc_Term2,vel_Term,delta_Term,dt);
    if (MODE==NVT) ABANH_update_cort(&gzi,&s,&s_vel,s_acc,predict_s,correct_s,dt);

    ABA_update(clt,crd,qrot,numclut,numatom);
    if (TERMMODE==ON) ABA_update_Term(crd,delta_Term,numatom);

    if (TERMMODE==ON)
      ABA_calcKineE_TermOn(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,s,s_vel,Q_NH,vel_Term,numclut,numatom,MODE,numclut+6);
    else
      ABA_calcKineE(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,s,s_vel,Q_NH,numclut,numatom,MODE);
    if (MODE==NVT)
      ABANH_calcKE(s,s_vel,Q_NH,KEobj,&PEv,&KEv);

    //    make_native_contact_list_aa(&num_NC,crd,numatom,criteria_NC,ncmap_aa,ON);
    //    Q_NC=((double)num_NC)/((double)e.num_natatt);

    //    Q_NC=count_native_contact(num_NC,crd,numatom,numres,indexncb,criteria_NC,ON);
    Q_NC=count_native_contact(num_NC,crd,numatom,numres,indexncb,criteria_NC,EXC);
    
    if (i%interval==0) {
      fprintf(outputfile,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i,pot,KE,KEv,PEv,pot+KE+PEv+KEv,s,T,Q_NC);
      //      GOLMAA_out_formated(outputfile2,e,KE,i,dt);
    }
    if (i%intervalnc==0) {
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crd[j*3+k];
      myncL_put_crd_AMBER(nc_id,l,crd_nc);
      ++l;
    }
  }
  fclose(outputfile);
  //  fclose(outputfile2);
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

  rstvelfile=efopen(rstvelfilename,"w");
  if (TERMMODE==ON) for (i=0;i<6;++i) fprintf(rstvelfile,"%10.8lf\n",vel_Term[i]);
  for (i=0;i<numclut;++i) fprintf(rstvelfile,"%10.8lf\n",qvel[i]);
  if (MODE==NVT) fprintf(rstvelfile,"%10.8lf %10.8lf",s,s_vel);
  fclose(rstvelfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfilename parmfilename outputfilename\n",progname);
}


