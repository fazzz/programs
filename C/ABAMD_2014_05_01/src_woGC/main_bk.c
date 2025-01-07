#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "netcdf.h"
#include <getopt.h>

//#include "ABA.h" // 2014-06-18
#include "ABAb.h"  // 2014-06-18

#include "PTL.h"
#include "FFL.h"
#include "EF.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

// 2014/04/23
//void check_maxabs_forceandTorque(int numatom,int numclut,double *frc,double *Q, 
//				 double *maxabs_force, int *index_atom, double *maxabs_Q, int *index_clust );
// 2014/04/23

// 2014/04/23
//// 2014/02/03
//void check_rmsd_forceandTorque(int numatom,int numclut,double *frc,double *Q, 
//			       double *rmsd_force,double *rmsd_Q );
//// 2014/02/03
// 2014/04/23

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d;
  int numatom,numclut,numstep=100000;
  int interval=1000,intervalout=1000,intervalnc=1000,intervalflag;
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

  double maxabs_force,maxabs_Q;
  int index_atom,index_clust;
  double maxabs_force_trj,maxabs_Q_trj;
  int index_atom_trj,index_clust_trj;

  //  // 2014/02/03                     // 2014-06-23
  //  double rmsd_force,rmsd_Q;         // 2014-06-23
  //  double ave_rmsd_force,ave_rmsd_Q; // 2014-06-23
  // 2014/02/03                         // 2014-06-23

  char *inputfilename,*reffilename,*outputfilename,*outputfilename2,*parmfilename,*clustfilename,*inputvelofilename;
  char *trjfilename;
  char *rstfilename="rstfile",*rstvelfilename="rstvelfile";
  FILE *inputfile,*reffile,*inputfilegoal,*outputfile,*outputfile2,*parmfile,*clustfile,*inputvelofile;
  FILE *velfile,*rstfile,*rstvelfile;

  //  char *check_forceandTorquefilename; // 2014-06-23
  //  FILE *check_forceandTorquefile;     // 2014-06-23

  // 2014/02/03                                  // 2014-06-23
  //  char *check_forceandTorquefilename_rmsd;   // 2014-06-23
  //  FILE *check_forceandTorquefile_rmsd;       // 2014-06-23
  // 2014/02/03                                  // 2014-06-23

  char *progname;

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
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"*+del14h@:t:k:c:t:o:s:i:j:k:v:?:{:}:",long_opt,&opt_idx))!=-1) {
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

  //  if (argc < /*6*/8/*2014/02/03*/) { // 2014-06-23
  if (argc < 6) {                        // 2014-06-23
    USAGE(progname);
    exit(1);
  }
  inputfilename     = *argv;
  clustfilename     = *++argv;
  parmfilename      = *++argv;
  outputfilename    = *++argv;
  outputfilename2   = *++argv;
  trjfilename       = *++argv;
  //  check_forceandTorquefilename = *++argv; // 2014-06-23
  // 2014/02/03                               // 2014-06-23
  //  check_forceandTorquefilename_rmsd = *++argv; // 2014-06-23
  
  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  clustfile=efopen(clustfilename,"r");
  //  clt=ABAp_clustscan(clustfile,&numclut); // 2014-06-20
  //  clt=(CLT *)ABAp_clustscan(clustfile,&numclut); // 2014-06-17
  //  clt=ABAp_clustscan(clustfile,&numclut); // 2014-06-19
  ABAp_clustscan(clustfile,&numclut,clt); // 2014-06-19
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
  ABANH_set_new(s,s_vel,gzi,predict_gzi,correct_gzi,predict_s,correct_s,tau,&tau2,&Q_NH,KEobj,dt);

  ABA_integ_set(q,qvel,predict,correct,numclut,dt);

  if (MODEV==ON) ABAs_restat_read_new(inputvelofilename,numclut,correct,correct_Term,correct_Term2,correct_gzi,MODE,TERMMODE);

  myncL_create_def_AMBER(trjfilename,numatom,&nc_id);
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");

  //check_forceandTorquefile=efopen(check_forceandTorquefilename,"w"); // 2014-06-23
  // 2014/02/03 // 2014-06-23
  //check_forceandTorquefile_rmsd=efopen(check_forceandTorquefilename_rmsd,"w"); // 2014-06-23
  // 2014/02/03 // 2014-06-23

  for (i=0;i<numstep;++i) {
    ABA_integ_pret(qrot,qvel,q,predict,correct,dt,numclut);
    if (TERMMODE==ON) ABA_integ_pret_Term(predict_Term,predict_Term2,correct_Term,correct_Term2,vel_Term,delta_Term,dt);

    if (MODE==NVT) ABANH_update_pret_new(&gzi,predict_gzi,correct_gzi,&s,&s_vel,predict_s,correct_s,dt);
    ABA_update(clt,crd,qrot,numclut,numatom);
    //    if (TERMMODE==ON) ABA_update_Term(crd,delta_Term,numatom);
    if (TERMMODE==ON) ABA_update_Term(crd,delta_Term,numatom,clt); // 2014-06-19

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

    //    check_maxabs_forceandTorque(numatom,numclut,frc,Q,&maxabs_force,&index_atom,&maxabs_Q,&index_clust);
    // 2014-06-23
    //    if (i==0) {                                // 2014-06-23
    //      maxabs_force_trj=maxabs_force;           // 2014-06-23
    //      index_atom_trj=index_atom;               // 2014-06-23
    //      maxabs_Q_trj=maxabs_Q;                   // 2014-06-23
    //      index_clust_trj=index_clust;             // 2014-06-23
    //    }                                          // 2014-06-23
    //    else {                                     // 2014-06-23
    //      if (maxabs_force_trj<maxabs_force) {     // 2014-06-23
    //	maxabs_force_trj=maxabs_force;               // 2014-06-23
    //	index_atom_trj=index_atom;                   // 2014-06-23
    //      }                                        // 2014-06-23
    //                                               // 2014-06-23
    //      if (maxabs_Q_trj<maxabs_Q) {             // 2014-06-23
    //	maxabs_Q_trj=maxabs_Q;                       // 2014-06-23
    //	index_clust_trj=index_clust;                 // 2014-06-23
    //      }                                        // 2014-06-23
    //    }                                          // 2014-06-23

    // 2014/02/03 // 2014-06-23
    //    check_rmsd_forceandTorque(numatom,numclut,frc,Q,&rmsd_force,&rmsd_Q); // 2014-06-23
    //    ave_rmsd_force=(ave_rmsd_force*i+rmsd_force)/(i+1);                   // 2014-06-23
    //    ave_rmsd_Q=(ave_rmsd_Q*i+rmsd_Q)/(i+1);                               // 2014-06-23
    // 2014/02/03

    ABA_calcKineE_TermOn(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,s,s_vel,Q_NH,vel_Term,numclut,numatom,MODE,numclut+6);
    T=KE/(DOF*k_B)*2.0;
    //    T=KE/(DOF*k_B)/*/2.0*/;
    if (TERMMODE==ON) 
      solverABA_TermOn_NH_new(qacc,qvel,clt,Q,frc,crd,numclut,numatom,q,gzi,&gzi_vel,s,&s_vel,tau2,acc_Term,acc_Term2,vel_Term,T,Tobj);
    else solverABA_NH_new(qacc,qvel,clt,Q,frc,crd,numclut,numatom,q,gzi,&gzi_vel,s,&s_vel,tau2,T,Tobj);

    ABA_integ_cort(qrot,qvel,q,qacc,predict,correct,dt,numclut);
    if (TERMMODE==ON)
      ABA_integ_cort_Term(predict_Term,predict_Term2,correct_Term,correct_Term2,acc_Term,acc_Term2,vel_Term,delta_Term,dt);
    if (MODE==NVT) ABANH_update_cort_new(&gzi,gzi_vel,&s,&s_vel,predict_gzi,correct_gzi,predict_s,correct_s,dt);

    ABA_update(clt,crd,qrot,numclut,numatom);
    //    if (TERMMODE==ON) ABA_update_Term(crd,delta_Term,numatom);
    if (TERMMODE==ON) ABA_update_Term(crd,delta_Term,numatom,clt); // 2014-06-19

    if (TERMMODE==ON)
      ABA_calcKineE_TermOn(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,s,s_vel,Q_NH,vel_Term,numclut,numatom,MODE,numclut+6);
    else
      ABA_calcKineE(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,s,s_vel,Q_NH,numclut,numatom,MODE);
    if (MODE==NVT)
      ABANH_calcKE_new(gzi,s,s_vel,Q_NH,KEobj,&PEv,&KEv);

    if (i%interval==0) 
      fprintf(outputfile,"%d %24.20e %24.20e %24.20e %24.20e %24.20e %24.20e %24.20e\n",i,pot,KE,KEv,PEv,pot+KE+PEv+KEv,s,T);

    if (i%intervalout==0) {
      //      ffL_out_formated(outputfile2,e,KE,KEv,PEv,T,i,dt);
      fprintf(outputfile2,"/***********************************************/\n");
      fprintf(outputfile2,"steps            = %d  \n",i);
      fprintf(outputfile2,"total time       = %10.3lf ps  \n" ,dt*(double)i);
      fprintf(outputfile2,"T_kelvin         = %24.20e K  \n",T);
      fprintf(outputfile2,"toal_energy      = %24.20e kcal/mol  \n",e.p_t+KE);
      fprintf(outputfile2,"toal_vertial_energy      = %24.20e kcal/mol  \n",e.p_t+KE+KEv+PEv);
      fprintf(outputfile2,"kinetic_energy   = %24.20e kcal/mol  \n",KE);
      fprintf(outputfile2,"kinetic_energy_vertial   = %24.20e kcal/mol  \n",KEv);
      fprintf(outputfile2,"potential_energy_real = %24.20e kcal/mol  \n",e.p_t);
      fprintf(outputfile2,"potential_energy_vertial   = %24.20e kcal/mol  \n",PEv);
      fprintf(outputfile2,"dihedral_energy  = %24.20e kcal/mol  \n",e.p_d_t);
      fprintf(outputfile2,"elect_energy     = %24.20e kcal/mol  \n",e.p_e_t);
      fprintf(outputfile2,"VDW_energy       = %24.20e kcal/mol  \n",e.p_LJ_t);
      fprintf(outputfile2,"1_4_elect_energy = %24.20e kcal/mol  \n",e.p_e_14_t);
      fprintf(outputfile2,"1_4_VDW_energy   = %24.20e kcal/mol  \n",e.p_LJ_14_t);

      //      fprintf(check_forceandTorquefile,"%4d %24.20e %4d %24.20e \n",index_atom+1,maxabs_force,index_clust+1,maxabs_Q); // 2014-06-23
      // 2014/02/03
      //      fprintf(check_forceandTorquefile_rmsd,"%24.20e %24.20e \n",rmsd_force,rmsd_Q); // 2014-06-23
      // 2014/02/03

    }

    if (i%intervalnc==0) {
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crd[j*3+k];
      myncL_put_crd_AMBER(nc_id,l,crd_nc);
      ++l;
    }
  }
  fclose(outputfile);
  fclose(outputfile2);

  //  fprintf(check_forceandTorquefile,"%4d %24.20e %4d %24.20e \n",index_atom_trj+1,maxabs_force_trj,index_clust_trj+1,maxabs_Q_trj); // 2014-06-23
  // 2014/02/03
  //  fprintf(check_forceandTorquefile_rmsd,"%24.20e %24.20e \n",ave_rmsd_force,ave_rmsd_Q); // 2014-06-23
  // 2014/02/03

  //  fclose(check_forceandTorquefile);      // 2014-06-23
  // 2014/02/03
  //  fclose(check_forceandTorquefile_rmsd); // 2014-06-23
  // 2014/02/03
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

  ABAs_restat_write_vel_new(rstvelfilename,numclut,correct,correct_Term,correct_Term2,correct_gzi,MODE,TERMMODE);
  //  fclose(outputfile);

  return 0;
}

//int USAGE(char *progname) {                                                             // 2-14-06-23
//  printf("USAGE:\n");                                                                   // 2-14-06-23
//  printf("[-h] help \n");                                                               // 2-14-06-23
//  printf("%s [-h] inputfilename clustfilename parmfilename outputfilename\n",progname); // 2-14-06-23
//}                                                                                       // 2-14-06-23

// 2-14-06-23
int USAGE(char *progname) {                                                             
  printf("USAGE:\n");                                                                   
  printf("[-f]\n");
  printf("[--nve]\n");
  printf("[--termon]\n");
  printf("[--dih]\n");
  printf("[--els]\n");
  printf("[--lj]\n");
  printf("[--e14]\n");
  printf("[--l14]\n");
  printf("[-h]\n");
  printf("[--kc]\n");
  printf("[--cof]\n");
  printf("[--temp]\n");
  printf("[--dt]\n");
  printf("[--omg]\n");
  printf("[--nums]\n");
  printf("[--int]\n");
  printf("[--intnc]\n");
  printf("[--intout]\n");
  printf("[--vel]\n");
  printf("[--tau]\n");
  printf("[--rst]\n");
  printf("[--rstv]\n");
  printf("[-h] help \n");                                                               
  printf("%s [-h] inputfilename clustfilename parmfilename outputfilename\n",progname); 
}                                                                                       
// 2014-06-23

// 2014-06-23
//void check_maxabs_forceandTorque(int numatom,int numclut,double *frc,double *Q, 
//				 double *maxabs_force, int *index_atom, double *maxabs_Q, int *index_clust ) {
//  int i,j;
//  double abs_force,abs_Q;
//
//  *maxabs_force=frc[0]*frc[0]+frc[1]*frc[1]+frc[2]*frc[2];
//  *maxabs_force=sqrt(*maxabs_force);
//  *index_atom=0;
//  for (i=1;i<numatom;++i) {
//    abs_force=frc[i*3]*frc[i*3]+frc[i*3+1]*frc[i*3+1]+frc[i*3+2]*frc[i*3+2];
//    abs_force=sqrt(abs_force);
//    if (abs_force>*maxabs_force) {
//      *maxabs_force=abs_force;
//      *index_atom=i;
//    }
//  }
//
//  *maxabs_Q=fabs(Q[0]);
//  *index_clust=0;
//  for (i=1;i<numclut;++i) {
//    abs_Q=fabs(Q[i]);
//    if (abs_Q>*maxabs_Q) {
//      *maxabs_Q=abs_Q;
//      *index_clust=i;
//    }
//  }
//
//}
// 2014-06-23

// 2014/02/03
//void check_rmsd_forceandTorque(int numatom,int numclut,double *frc,double *Q, 
//			       double *rmsd_force,double *rmsd_Q ) {
//  int i,j;
//  double ave_force[3],ave_Q;
//
//  for (i=1;i<3;++i) ave_force[i]=0.0;
//
//  for (i=1;i<numatom;++i)
//    for (j=1;j<3;++j)  
//      ave_force[j]+=frc[i*3+j];
//
//  for (i=1;i<3;++i)  ave_force[i]/=numatom;
//
//  *rmsd_force=0.0;
//  for (i=1;i<numatom;++i) {
//    for (j=1;j<3;++j) {
//      *rmsd_force+=(frc[i*3+j]-ave_force[j])*(frc[i*3+j]-ave_force[j]);
//    }
//  }
//  *rmsd_force=*rmsd_force/numatom;
//  *rmsd_force=sqrt(*rmsd_force);
//
//  ave_Q=0.0;
//  for (i=0;i<numclut;++i) ave_Q+=Q[i];
//  ave_Q/=numclut;
//
//  *rmsd_Q=0.0;
//  for (i=0;i<numclut;++i) {
//    *rmsd_Q+=(Q[i]-ave_Q)*(Q[i]-ave_Q);
//  }
//  *rmsd_Q=*rmsd_Q/numclut;
//  *rmsd_Q=sqrt(*rmsd_Q);
//
//}
// 2014/02/03
