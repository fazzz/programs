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
  int numatom,numclut,numstep=100000;
  int interval=1000,intervalout=1000,intervalnc=1000,intervalflag;
  double dt=0.001;
  CLTb *clt;
  double *Q,*frc,pot;
  double *q;
  double *qacc,*qvel,*qrot;
  double *predict,*correct;
  double KE,KEv,PEv;
  double dummy;

  int MODE=NVT,MODEV=OFF,TERMMODE=OFF,TUMODE=OFF;
  double Tobj=300,KEobj;
  double k_B=1.98723e-3;

  int *atom_tune_pairs;
  double /**tune_val_by_period*/*tune_val_V_n,*tune_val_n_phase;
  int numtune;

  /***************************/
  /* int *atom_tune_14pairs; */
  /* double *tune_14val;     */
  /* int num14tune;	     */
  /***************************/

  int *atom_tune_14pairs_es,*atom_tune_14pairs_LJ;
  double *tune_14val_es,*tune_14val_LJ;
  int num14tune_es,num14tune_LJ;

  double *p_d,p_d_t;

  double s=1.0,s_vel=0.0,s_acc,gzi=0.0,gzi_vel,predict_gzi[5],correct_gzi[5],predict_s[5],correct_s[5];
  double Q_NH,tau=0.1,tau2;
  double T;
  int DOF;

  int dfalg=ON,efalg=ON,LJfalg=ON,e14falg=ON,LJ14falg=ON;
  double coff=1.0;
  int *numclutparent,*terminal,*origin;

  double *delta_Term,*vel_Term,*acc_Term,*acc_Term2,**predict_Term,**predict_Term2,**correct_Term,**correct_Term2;

  double /***crd_nc,*/*energy;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_AMBER nc_id;

  double *crd,*mass;
  struct potential e;
  struct force f;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*reffilename,*outputfilename,*outputfilename2,*parmfilename,*clustfilename,*inputvelofilename;
  char *trjfilename;
  char *rstfilename="rstfile",*rstvelfilename="rstvelfile",*TUMEDIHEDfilename;
  FILE *inputfile,*reffile,*inputfilegoal,*outputfile,*outputfile2,*parmfile,*clustfile,*inputvelofile;
  FILE *velfile,*rstfile,*rstvelfile,*TUMEDIHEDfile;

  char *progname;

  double pi;

  int opt_idx=1;

  struct option long_opt[] = {
    {"f",0,NULL,'f'},
    {"nve",0,NULL,'*'},
    {"termon",0,NULL,'+'},
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
    {"intnc",1,NULL,'j'},
    {"intout",1,NULL,'k'},
    {"vel",1,NULL,'v'},
    {"tau",1,NULL,'?'},
    {"rst",1,NULL,'{'},
    {"rstv",1,NULL,'}'},
    {"tune",1,NULL,'T'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"*+del14h@:t:k:c:t:o:s:i:j:k:v:?:{:}:T:",long_opt,&opt_idx))!=-1) {
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
    case 'e':
      efalg=OFF;
      break;
    case 'l':
      LJfalg=OFF;
      break;
    case '1':
      e14falg=OFF;
      break;
    case '4':
      LJ14falg=OFF;
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
    case 'c':
      coff=atof(optarg);
      break;
    case 'v':
      inputvelofilename=optarg;
      MODEV=ON;
      break;
    case '?':
      tau=atof(optarg);
      break;
    case '{':
      rstfilename=optarg;
      break;
    case '}':
      rstvelfilename=optarg;
      break;
    case 'T':
      TUMODE=ON;
      TUMEDIHEDfilename=optarg;
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  pi=acos(-1.0);

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
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

  if (TUMODE==ON) {
    TUMEDIHEDfile=efopen(TUMEDIHEDfilename,"r");
    fscanf(TUMEDIHEDfile,"%d",&numtune);
    atom_tune_pairs=(int *)gcemalloc(sizeof(int)*4*numtune);
    tune_val_V_n=(double *)gcemalloc(sizeof(double)*numtune);
    tune_val_n_phase=(double *)gcemalloc(sizeof(double)*numtune);
    for (i=0;i<numtune;++i) {
      for (j=0;j<4;++j)	fscanf(TUMEDIHEDfile,"%d",&atom_tune_pairs[i*4+j]);
      for (j=0;j<4;++j) atom_tune_pairs[i*4+j]-=1;
      fscanf(TUMEDIHEDfile,"%lf",&tune_val_V_n[i]);
      fscanf(TUMEDIHEDfile,"%lf",&tune_val_n_phase[i]);
    }
    fscanf(TUMEDIHEDfile,"%d",&num14tune_es);
    for (i=0;i<num14tune_es;++i) {
      atom_tune_14pairs_es=(int *)gcemalloc(sizeof(int)*2*num14tune_es);
      tune_14val_es=(double *)gcemalloc(sizeof(double)*num14tune_es);
      for (j=0;j<2;++j)	fscanf(TUMEDIHEDfile,"%d",&atom_tune_14pairs_es[i*2+j]);
      for (j=0;j<2;++j) atom_tune_14pairs_es[i*2+j]-=1;
      fscanf(TUMEDIHEDfile,"%lf",&tune_14val_es[i]);
    }
    fscanf(TUMEDIHEDfile,"%d",&num14tune_LJ);
    for (i=0;i<num14tune_LJ;++i) {
      atom_tune_14pairs_LJ=(int *)gcemalloc(sizeof(int)*2*num14tune_LJ);
      tune_14val_LJ=(double *)gcemalloc(sizeof(double)*num14tune_LJ);
      for (j=0;j<2;++j)	fscanf(TUMEDIHEDfile,"%d",&atom_tune_14pairs_LJ[i*2+j]);
      for (j=0;j<2;++j) atom_tune_14pairs_LJ[i*2+j]-=1;
      fscanf(TUMEDIHEDfile,"%lf",&tune_14val_LJ[i]);
    }
    fclose(TUMEDIHEDfile);
    p_d=(double *)gcemalloc(sizeof(double)*(AP.NPHIH+AP.MPHIA));
  }

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

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

  q=(double *)gcemalloc(sizeof(double)*numclut);
  qacc=(double *)gcemalloc(sizeof(double)*numclut);
  qvel=(double *)gcemalloc(sizeof(double)*numclut);
  qrot=(double *)gcemalloc(sizeof(double)*numclut);
  predict=(double *)gcemalloc(sizeof(double)*numclut*6/**6*/);
  correct=(double *)gcemalloc(sizeof(double)*numclut*6/**6*/);

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
  for (i=0;i<numclut;++i) qvel[i]=0.00/*1.00*/;

  DOF=(numclut-1);
  if (TERMMODE==ON) DOF+=6;
  KEobj=0.5*DOF*k_B*Tobj;
  //  KEobj=/*0.5**/DOF*k_B*Tobj;
  ABAbNH_set_new(s,s_vel,gzi,predict_gzi,correct_gzi,predict_s,correct_s,tau,&tau2,&Q_NH,KEobj,dt);

  ABAb_integ_set(q,qvel,predict,correct,numclut,dt);

  if (MODEV==ON) ABAbs_restat_read_new(inputvelofilename,numclut,correct,correct_Term,correct_Term2,correct_gzi,MODE,TERMMODE);

  myncL_create_def_AMBER(trjfilename,numatom,&nc_id);
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");

  for (i=0;i<numstep;++i) {
    ABAb_integ_pret(qrot,qvel,q,predict,correct,dt,numclut);
    if (TERMMODE==ON) ABAb_integ_pret_Term(predict_Term,predict_Term2,correct_Term,correct_Term2,vel_Term,delta_Term,dt);

    if (MODE==NVT) ABAbNH_update_pret_new(&gzi,predict_gzi,correct_gzi,&s,&s_vel,predict_s,correct_s,dt);
    ABAb_update(clt,crd,qrot,numclut,numatom);
    //    if (TERMMODE==ON) ABAb_update_Term(crd,delta_Term,numatom);

    for (j=0;j<numclut;++j) Q[j]=0.0;
    intervalflag=i%interval;
    if (dfalg==ON) {
      if (/*TERMMODE==ON*/TUMODE==ON)
	//	p_d_t=ffL_calcTorque_wtune(Q,p_d,crd,numclut,numclutparent,terminal,origin,atom_tune_pairs,tune_val_by_period,numtune);
	p_d_t=ffL_calcTorque_wtuneb(Q,p_d,crd,numclut,numclutparent,terminal,origin,atom_tune_pairs,tune_val_V_n,tune_val_n_phase,numtune,pi);
      else
	ffL_calcTorque(Q,crd,numclut,numclutparent,terminal,origin);
    }
    if (TUMODE==ON) {
      //      ffL_calcffandforce_w14tune(crd,numatom,&e,&f,atom_tune_14pairs,tune_14val,num14tune);
      ffL_calcffandforce_w14tuneb(crd,numatom,&e,&f,atom_tune_14pairs_es,atom_tune_14pairs_LJ,tune_14val_es,tune_14val_LJ,num14tune_es,num14tune_LJ);
    }
    else
      ffL_calcffandforce(crd,numatom,&e,&f);
    pot=0.0;
    if (dfalg==ON) {
      if (TERMMODE==ON)
	pot=p_d_t;
      else
	pot=e.p_d_t;
    }
    for (j=0;j<numatom;++j) for (k=0;k<3;++k) frc[j*3+k]=0.0;
    for (j=0;j<numatom;++j) {
      if (efalg   == ON) {
    	for (k=0;k<3;++k) frc[j*3+k]+=/*coff**/f.f_e[j*3+k];
      }
      if (LJfalg  == ON) {
    	for (k=0;k<3;++k) frc[j*3+k]+=/*coff**/f.f_LJ[j*3+k];
      }
      if (e14falg == ON) {
    	for (k=0;k<3;++k) frc[j*3+k]+=/*coff**/f.f_e_14[j*3+k];
      }
      if (LJ14falg== ON) {
    	for (k=0;k<3;++k) frc[j*3+k]+=/*coff**/f.f_LJ_14[j*3+k];
      }
    }
    if (efalg   == ON)  pot+=0.5*e.p_e_t;
    if (LJfalg  == ON)  pot+=0.5*e.p_LJ_t;
    if (e14falg == ON)  pot+=0.5*e.p_e_14_t;
    if (LJ14falg== ON)  pot+=0.5*e.p_LJ_14_t;

    ABAb_calcKineE_TermOn(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,s,s_vel,Q_NH,vel_Term,numclut,numatom,MODE,numclut+6);
    T=KE/(DOF*k_B)*2.0;
    //    T=KE/(DOF*k_B)/*/2.0*/;
    if (TERMMODE==ON) 
      solverABAb_TermOn_NH_new(qacc,qvel,clt,Q,frc,crd,numclut,numatom,q,gzi,&gzi_vel,s,&s_vel,tau2,acc_Term,acc_Term2,vel_Term,T,Tobj);
    else solverABAb_NH_new(qacc,qvel,clt,Q,frc,crd,numclut,numatom,q,gzi,&gzi_vel,s,&s_vel,tau2,T,Tobj);

    ABAb_integ_cort(qrot,qvel,q,qacc,predict,correct,dt,numclut);
    if (TERMMODE==ON)
      ABAb_integ_cort_Term(predict_Term,predict_Term2,correct_Term,correct_Term2,acc_Term,acc_Term2,vel_Term,delta_Term,dt);
    if (MODE==NVT) ABAbNH_update_cort_new(&gzi,gzi_vel,&s,&s_vel,predict_gzi,correct_gzi,predict_s,correct_s,dt);

    ABAb_update(clt,crd,qrot,numclut,numatom);
    //    if (TERMMODE==ON) ABAb_update_Term(crd,delta_Term,numatom);

    if (TERMMODE==ON)
      ABAb_calcKineE_TermOn(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,s,s_vel,Q_NH,vel_Term,numclut,numatom,MODE,numclut+6);
    else
      ABAb_calcKineE(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,s,s_vel,Q_NH,numclut,numatom,MODE);
    if (MODE==NVT)
      ABAbNH_calcKE_new(gzi,s,s_vel,Q_NH,KEobj,&PEv,&KEv);

    if (i%interval==0) 
      fprintf(outputfile,"%d %e %e %e %e %e %e %e\n",i,pot,KE,KEv,PEv,pot+KE+PEv+KEv,s,T);

    if (i%intervalout==0) {
      //      ffL_out_formated(outputfile2,e,KE,KEv,PEv,T,i,dt);
      fprintf(outputfile2,"/***********************************************/\n");
      fprintf(outputfile2,"steps            = %d  \n",i);
      fprintf(outputfile2,"total time       = %10.3lf ps  \n" ,dt*(double)i);
      fprintf(outputfile2,"T_kelvin         = %e K  \n",T);
      fprintf(outputfile2,"toal_energy      = %e kcal/mol  \n",e.p_t+KE);
      fprintf(outputfile2,"toal_vertial_energy      = %e kcal/mol  \n",e.p_t+KE+KEv+PEv);
      fprintf(outputfile2,"kinetic_energy   = %e kcal/mol  \n",KE);
      fprintf(outputfile2,"kinetic_energy_vertial   = %e kcal/mol  \n",KEv);
      fprintf(outputfile2,"potential_energy_real = %e kcal/mol  \n",e.p_t);
      fprintf(outputfile2,"potential_energy_vertial   = %e kcal/mol  \n",PEv);
      fprintf(outputfile2,"dihedral_energy  = %e kcal/mol  \n",e.p_d_t);
      fprintf(outputfile2,"elect_energy     = %e kcal/mol  \n",e.p_e_t);
      fprintf(outputfile2,"VDW_energy       = %e kcal/mol  \n",e.p_LJ_t);
      fprintf(outputfile2,"1_4_elect_energy = %e kcal/mol  \n",e.p_e_14_t);
      fprintf(outputfile2,"1_4_VDW_energy   = %e kcal/mol  \n",e.p_LJ_14_t);
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

  ABAbs_restat_write_vel_new(rstvelfilename,numclut,correct,correct_Term,correct_Term2,correct_gzi,MODE,TERMMODE);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfilename parmfilename outputfilename\n",progname);
}


