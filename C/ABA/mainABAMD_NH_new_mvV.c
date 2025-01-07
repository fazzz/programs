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

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d,t;
  int numatom,numclut,numstep=100000;
  int interval=1000,intervalout=1000,intervalnc=1000,intervalflag,numite=2;
  double dt=0.001;
  CLT *clt;
  double *Q,*frc,pot;
  double *q;
  double *qacc,*qvel,*qrot;
  double *qvel_b1,*qvel_b2;
  double *qvel_b1_temp,*qvel_b2_temp;
  double KE,KE_1_2,KEv,PEv,KEv_1_2,PEv_1_2;
  double dummy;

  int MODE=NVT,MODEV=OFF,TERMMODE=OFF;
  double Tobj=300,KEobj;
  double k_B=1.98723e-3;

  double s=1.0,s_vel=0.0,s_acc,zeta=0.0,zeta_vel;
  double Q_NH,tau=0.1,tau2;
  double T,T_1_2;
  int DOF;

  int dfalg=ON,efalg=ON,LJfalg=ON,e14falg=ON,LJ14falg=ON;
  double coff=1.0;
  int *numclutparent,*terminal,*origin;

  double *delta_Term,*vel_Term,*acc_Term,*acc_Term2;
  double *vel_Term_b1,*vel_Term_b2;
  double *vel_Term_b1_temp,*vel_Term_b2_temp;

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
    {"els",0,NULL,'e'},
    {"lj",0,NULL,'l'},
    {"e14",0,NULL,'1'},
    {"l14",0,NULL,'4'},
    {"h",0,NULL,'h'},
    //    {"kc",1,NULL,'k'},
    {"cof",1,NULL,'c'},
    {"temp",1,NULL,'t'},
    {"dt",1,NULL,'@'},
    {"omg",1,NULL,'o'},
    {"nums",1,NULL,'s'},
    {"int",1,NULL,'i'},
    {"intnc",1,NULL,'j'},
    {"intout",1,NULL,'k'},
    {"ite",1,NULL,'x'},
    {"vel",1,NULL,'v'},
    {"tau",1,NULL,'?'},
    {"rst",1,NULL,'{'},
    {"rstv",1,NULL,'}'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"*+del14h@:t:c:t:o:s:i:j:k:x:v:?:{:}:",long_opt,&opt_idx))!=-1) {
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
    case 'x':
      numite=atoi(optarg);
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

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

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
  qrot=(double *)gcemalloc(sizeof(double)*numclut);
  qvel_b1=(double *)gcemalloc(sizeof(double)*numclut*6);
  qvel_b2=(double *)gcemalloc(sizeof(double)*numclut*6);
  qvel_b1_temp=(double *)gcemalloc(sizeof(double)*numclut*6);

  delta_Term=(double *)gcemalloc(sizeof(double)*6);
  vel_Term=(double *)gcemalloc(sizeof(double)*6);
  acc_Term=(double *)gcemalloc(sizeof(double)*6);
  acc_Term2=(double *)gcemalloc(sizeof(double)*6);
  vel_Term_b1=(double *)gcemalloc(sizeof(double)*6);
  vel_Term_b2=(double *)gcemalloc(sizeof(double)*6);
  vel_Term_b1_temp=(double *)gcemalloc(sizeof(double)*6);

  for (i=0;i<6;++i) {
    vel_Term[i]=0.0;
    vel_Term_b1[i]=0.0;
    vel_Term_b2[i]=0.0;
  }
  for (i=0;i<numclut;++i) {
    qvel[i]=0.00/*1.00*/;
    qvel_b1[i]=0.00;
    qvel_b2[i]=0.00;
  }

  DOF=(numclut-1);
  if (TERMMODE==ON) DOF+=6;
  KEobj=0.5*DOF*k_B*Tobj;

  ABANH_set_new_mvV(s,s_vel,zeta,tau,&tau2,&Q_NH,KEobj,dt);

  if (MODEV==ON) 
    ABAs_restat_read_new_mvV(inputvelofilename,numclut,qvel_b1,qvel_b2,vel_Term_b1,vel_Term_b2,zeta,MODE,TERMMODE);

  myncL_create_def_AMBER(trjfilename,numatom,&nc_id);
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");

  for (i=0;i<numstep;++i) {
    for (j=0;j<numclut;++j) qvel[j]=1.5*qvel_b1[j]-0.5*qvel_b2[j];
    if (TERMMODE==ON) for (j=0;j<6;++j) vel_Term[j]=1.5*vel_Term_b1[j]-0.5*vel_Term_b2[j];

    
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
    
    for (t=0;t<numite;++t) {
      if (TERMMODE==ON) 
	solverABA_TermOn_NH_new_mvV(qacc,qvel,clt,Q,frc,crd,numclut,numatom,q,zeta,acc_Term,acc_Term2,vel_Term);
      else solverABA_NH_new_mvV(qacc,qvel,clt,Q,frc,crd,numclut,numatom,q,zeta);
      
      for (j=0;j<numclut;++j) {
	/*********************************************/
        /* qvel[j]=qvel_b2[j]+0.5*dt*qacc[j];	     */
	/* qvel_b2[j]=qvel_b1[j];		     */
	/* qvel_b1[j]=qvel_b1[j]+dt*qacc[j];	     */
        /*********************************************/
	//      qrot[j]=dt*qvel[j]+0.5*dt*dt*qacc[j];
	/*********************************/
        /* qrot[j]=dt*qvel_b1[j];	 */
        /*********************************/
	qvel_b1_temp[j]=qvel_b1[j]+dt*qacc[j];
	qvel[j]=0.5*qvel_b1[j]+0.5*qvel_b1_temp[j];
	//	qrot[j]=dt*qvel_b1_temp[j];
	//	qrot[j]=dt*qvel[j]+0.5*dt*dt*qacc[j];
      }
    
      if (TERMMODE==ON) {
	for (j=0;j<6;++j) {
	  /********************************************************************/
          /* vel_Term[j]=vel_Term_b2[j]+0.5*dt*acc_Term[j];		      */
	  /* vel_Term_b2[j]=vel_Term_b1[j];				      */
	  /* vel_Term_b1[j]=vel_Term_b1[j]+dt*acc_Term[j];		      */
	  /* //	delta_Term[j]=dt*vel_Term[j]+0.5*dt*dt*acc_Term[j];	      */
	  /* delta_Term[j]=dt*vel_Term_b1[j];				      */
          /********************************************************************/
	  vel_Term_b1_temp[j]=vel_Term_b1[j]+dt*acc_Term[j];
	  vel_Term[j]=0.5*vel_Term_b1[j]+0.5*vel_Term_b1_temp[j];
	  //	  delta_Term[j]=dt*vel_Term_b1_temp[j];
	  //	  delta_Term[j]=dt*vel_Term[j]+0.5*dt*dt*acc_Term[j];
	}
      }
    }

    for (j=0;j<numclut;++j) {
      qrot[j]=/*dt*qvel_b1_temp[j]*/dt*qvel[j]+0.5*dt*dt*qacc[j];
      qvel_b2[j]=qvel_b1[j];
      qvel_b1[j]=qvel_b1_temp[j];
    }
    if (TERMMODE==ON) {
      for (j=0;j<6;++j) {
	delta_Term[j]=/*dt*vel_Term_b1_temp[j]*/dt*vel_Term[j]+0.5*dt*dt*acc_Term[j];
	vel_Term_b2[j]=vel_Term_b1[j];
	vel_Term_b1[j]=vel_Term_b1_temp[j];
      }
    }
      
    if (MODE==NVT) { 
      ABA_calcKineE_TermOn(&KE_1_2,&KEv_1_2,&PEv_1_2,KEobj,clt,crd,qvel_b1,s,s_vel,Q_NH,vel_Term_b1,numclut,numatom,MODE,numclut+6);
      T_1_2=KE_1_2/(DOF*k_B)*2.0;
      zeta=zeta+dt*1.0/tau2*(T_1_2/Tobj-1.0);
      //////////////////////////////////////
      s_vel=zeta*s;
      s=s+dt*s_vel;
      //////////////////////////////////////
      }
    
    ABA_update(clt,crd,qrot,numclut,numatom);
    //    if (TERMMODE==ON) ABA_update_Term(crd,delta_Term,numatom);
    //  }

    if (TERMMODE==ON)
      ABA_calcKineE_TermOn(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,s,s_vel,Q_NH,vel_Term,numclut,numatom,MODE,numclut+6);
    else
      ABA_calcKineE(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,s,s_vel,Q_NH,numclut,numatom,MODE);
    if (MODE==NVT)
      ABANH_calcKE_new(zeta,s,s_vel,Q_NH,KEobj,&PEv,&KEv);
    T=KE/(DOF*k_B)*2.0;

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

  ABAs_restat_write_vel_new_mvV(rstvelfilename,numclut,qvel_b1,qvel_b2,vel_Term_b1,vel_Term_b2,zeta,MODE,TERMMODE);
  fclose(outputfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfilename parmfilename outputfilename\n",progname);
}


