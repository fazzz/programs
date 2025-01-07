#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "ABAb.h"
#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"
#include "GOLMAA_MB_PROTEINS2008.h"

#include "PTL.h"
#include "FFL.h"
#include "EF.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,dummyd;
  int joinflag;
  int numatom,numheavyatom,numres,numclut,numstep=10000,interval=100,interval2=100,interval3=100;
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
  double UNITT=418.4070;
  int forceflag=ON;

  double avePE=0.0,varPE=0.0,aveKE=0.0,varKE=0.0,aveT=0.0,varT=0.0;

  //  double *zeta,*zeta_vel,*zeta_acc,**predict_zeta,**correct_zeta;
  //  double *Q_NH,tau=0.1,tau2;
  //  double T;
  //  int DOF,M=4;
  double s=1.0,s_vel=0.0,s_acc,gzi=0.0,gzi_vel,predict_gzi[5],correct_gzi[5],predict_s[5],correct_s[5];
  double Q_NH,tau=0.1,tau2;
  double T,T_1_2;
  int DOF;

  double ep=0.3;
  int NCmode=3,nibnum=3,criteria=6.5;

  int *numclutparent,*terminal,*origin;
  double *delta_Term,*vel_Term,*acc_Term,*acc_Term2,**predict_Term,**predict_Term2,**correct_Term,**correct_Term2;

  double /***crd_nc,*/*energy;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_AMBER nc_id;

  double *crd,*refcrd1,*refcrd2,*mass;
  double d=1.0,de=1.0,d2;
  int KBTUNIT=ON;

  double summass,COM[3];

  struct potential e;
  struct force f;
  struct potential_GOLMAA_PROTEINS2008 e_GOLM_test;
  struct potential_GOLMAA_MB_PROTEINS2008 e_GOLM;
  double p_t=0.0,E_t;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*reffilename1,*reffilename2,*inputvelofilename,*parmfilename,*clustfilename;
  char *trjfilename,*outputfilename,*outputfilename2,*rstfilename="rstcrd",*rstvelfilename="rstvel";

  char *logfilename="MD_pep_NH_MP1996_GOLMAA_PROTEINS2008.log";

  FILE *inputfile,*reffile1,*reffile2,*velfile,*parmfile,*clustfile;
  FILE *outputfile,*outputfile2,*rstfile,*rstvelfile;

  FILE *logfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"nve",0,NULL,'*'},
    {"f",0,NULL,'f'},
    {"termon",0,NULL,'+'},
    {"vel",1,NULL,'v'},
    {"tau",1,NULL,'o'},
    {"nchain",1,NULL,'M'},
    {"ep",1,NULL,'e'},
    {"nums",1,NULL,'@'},
    {"temp",1,NULL,'t'},
    {"int",1,NULL,'l'},
    {"int1",1,NULL,'i'},
    {"int2",1,NULL,'j'},
    {"int3",1,NULL,'k'},
    {"rst",1,NULL,'{'},
    {"rstvel",1,NULL,'}'},
    {"dt",1,NULL,'x'},
    {"cutoff",1,NULL,'c'},
    {"de",1,NULL,'d'},
    {"dh",1,NULL,'2'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"*fh+ovM:e:@:t:l:i:j:k:{:}:x:c:d:2:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'f':
      forceflag=OFF;
      break;
    case '*':
      MODE=NVE;
      break;
    case '+':
      TERMMODE=ON;
      break;
    case 'o':
      tau=atof(optarg);
      break;
    /*********************/
    /* case 'M':	 */
    /*   M=atoi(optarg); */
    /*   break;		 */
    /*********************/
    case '@':
      numstep=atoi(optarg);
      break;
    case 'e':
      ep=atof(optarg);
      break;
    case 't':
      Tobj=atof(optarg);
      break;
    case 'x':
      dt=atof(optarg);
      break;
    case 'l':
      interval=atoi(optarg);
      interval2=interval;
      interval3=interval;
      break;
    case 'i':
      interval=atoi(optarg);
      break;
    case 'j':
      interval2=atoi(optarg);
      break;
    case 'k':
      interval3=atoi(optarg);
      break;
    case 'v':
      inputvelofilename=optarg;
      MODEV=ON;
      break;
    case '{':
      rstfilename=optarg;
      break;
    case '}':
      rstvelfilename=optarg;
      break;
    case 'c':
      criteria=atof(optarg);
      break;
    case 'd':
      de=atof(optarg);
      break;
    case '2':
      d=atof(optarg);
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

  if (argc < 8) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  reffilename1 = *++argv;
  reffilename2 = *++argv;
  clustfilename  = *++argv;
  parmfilename   = *++argv;
  outputfilename = *++argv;
  outputfilename2= *++argv;
  trjfilename    = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;

  j=0;
  for (i=0;i<numatom;++i) {
    if (strncmp(AP.IGRAPH[i],"H",1)==0) {
      ++j;
    }
  }
  numheavyatom=numatom-j;
  numres=AP.NRES;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];
  
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd1=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd2=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&dummyd);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  reffile1=efopen(reffilename1,"r");
  getline(&line,&len,reffile1);
  fscanf(reffile1,"%d",&dummyd);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(reffile1,"%lf",&refcrd1[i*3+j]);
  fclose(reffile1);

  reffile2=efopen(reffilename2,"r");
  getline(&line,&len,reffile2);
  fscanf(reffile2,"%d",&dummyd);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(reffile2,"%lf",&refcrd2[i*3+j]);
  fclose(reffile2);

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

  clt[0].join=0;
  for(i=1; i<numclut; ++i) {
    ABAb_setJoin(clt,i);
  }
  ABAbs_local_reference(clt,numclut,numatom,crd);
  ABAbs_trans_Matrix(clt,numclut,numatom,crd);
  ABAbs_inertia_matrix(clt,numclut,numatom,crd,mass);

  Q=(double *)gcemalloc(sizeof(double)*numclut);
  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

  summass=0.0;
  for (i=0;i<numatom;++i) summass+=mass[i];

  ffL_set_calcffandforce(&e,&f);
  GOLMAA_MB_PROTEINS2008_ff_calcff_set(&e_GOLM,refcrd1,refcrd2,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);


  /**************************************************************************************************************************/
  /* GOLMAA_PROTEINS2008_ff_set_calcff(&e_GOLM_test,refcrd1,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria); */
  /* GOLMAA_PROTEINS2008_ff_calcff_wobaimp(crd,numatom,&e_GOLM_test);							    */
  /**************************************************************************************************************************/

  q=(double *)gcemalloc(sizeof(double)*numclut);
  qacc=(double *)gcemalloc(sizeof(double)*numclut);
  qvel=(double *)gcemalloc(sizeof(double)*numclut);
  qrot=(double *)gcemalloc(sizeof(double)*numclut);
  predict=(double *)gcemalloc(sizeof(double)*numclut*6/**6*/);
  correct=(double *)gcemalloc(sizeof(double)*numclut*6/**6*/);

  /************************************************************/
  /* zeta=(double *)gcemalloc(sizeof(double)*M);	      */
  /* zeta_vel=(double *)gcemalloc(sizeof(double)*M);	      */
  /* zeta_acc=(double *)gcemalloc(sizeof(double)*M);	      */
  /* predict_zeta=(double **)gcemalloc(sizeof(double *)*M);   */
  /* correct_zeta=(double **)gcemalloc(sizeof(double *)*M);   */
  /* for (i=0;i<M;++i) {				      */
  /*   predict_zeta[i]=(double *)gcemalloc(sizeof(double)*6); */
  /*   correct_zeta[i]=(double *)gcemalloc(sizeof(double)*6); */
  /* }							      */
  /* Q_NH=(double *)gcemalloc(sizeof(double)*M);	      */
  /************************************************************/

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
  // debug
  for (i=0;i<numclut;++i) qvel[i]=/*0.10*//*0.01*/0.00;
  // debug

  DOF=(numclut-1);
  if (TERMMODE==ON) DOF+=6;
  KEobj=0.5*DOF*k_B*Tobj;
  KBT=k_B*Tobj;

  if (KBTUNIT==ON) {
    de=de*KBT;
    d=d*KBT;
  }
  
  d2=d*d;
  GOLMAA_MB_PROTEINS2008_ff_calcff_wobaimp(crd,numatom,de,d2,&e_GOLM);


  ABAbNH_set_new(s,s_vel,gzi,predict_gzi,correct_gzi,predict_s,correct_s,tau,&tau2,&Q_NH,KEobj,dt);
  //  ABAbNH_set_new(s,s_vel,gzi,predict_gzi,correct_gzi,predict_s,correct_s,tau,&tau2,&Q_NH,KEobj,dt);

  ABAb_integ_set(q,qvel,predict,correct,numclut,dt);

  if (MODEV==OFF) 
    ;
    //    KE=ABAb_Generate_inivelo(clt,crd,qvel,vel_Term,numclut,numatom,TERMMODE,KBT*UNITT);
  else
    ABAbs_restat_read_new(inputvelofilename,numclut,correct,correct_Term,correct_Term2,correct_gzi,MODE,TERMMODE);
  //    ABAbs_restat_read_chain(inputvelofilename,numclut,correct,correct_Term,correct_Term2,correct_zeta,M,MODE,TERMMODE);
  
  myncL_create_def_AMBER(trjfilename,numatom,&nc_id);
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");

  for (i=0;i<numstep;++i) {
    ABAb_integ_pret(qrot,qvel,q,predict,correct,dt,numclut);
    if (TERMMODE==ON) 
      ABAb_integ_pret_Term(predict_Term,predict_Term2,correct_Term,correct_Term2,vel_Term,delta_Term,dt);
    if (MODE==NVT) 
      ABAbNH_update_pret_new(&gzi,predict_gzi,correct_gzi,&s,&s_vel,predict_s,correct_s,dt);
    //      ABAbNH_chain_update_pret(zeta,zeta_vel,predict_zeta,correct_zeta,M,dt);

    ABAb_update(clt,crd,qrot,numclut,numatom);

    for (j=0;j<numclut;++j) Q[j]=0.0;
    for (j=0;j<numatom;++j) for (k=0;k<3;++k) frc[j*3+k]=0.0;

    GOLMAA_MB_PROTEINS2008_ff_calcff_wobaimp(crd,numatom,de,d2,&e_GOLM);
    for (j=0;j<numatom;++j)
      for (k=0;k<3;++k)
	frc[j*3+k]=e_GOLM.f_MB[j][k];

    ABAb_calcKineE_TermOn(&KE,&KEv,&PEv,KEobj,clt,crd,qvel,s,s_vel,Q_NH,vel_Term,numclut,numatom,MODE,numclut+6);
    T=KE/(DOF*k_B)*2.0;

    if (TERMMODE==ON) 
      solverABAb_TermOn_NH_new(qacc,qvel,clt,Q,frc,crd,numclut,numatom,q,gzi,&gzi_vel,s,&s_vel,tau2,acc_Term,acc_Term2,vel_Term,T,Tobj);
    else 
      solverABAb_NH_new(qacc,qvel,clt,Q,frc,crd,numclut,numatom,q,gzi,&gzi_vel,s,&s_vel,tau2,T,Tobj);

    ABAb_integ_cort(qrot,qvel,q,qacc,predict,correct,dt,numclut);
    if (TERMMODE==ON)
      ABAb_integ_cort_Term(predict_Term,predict_Term2,correct_Term,correct_Term2,acc_Term,acc_Term2,vel_Term,delta_Term,dt);
    if (MODE==NVT) ABAbNH_update_cort_new(&gzi,gzi_vel,&s,&s_vel,predict_gzi,correct_gzi,predict_s,correct_s,dt);
    //    ABAbNH_chain_update_cort(zeta,zeta_vel,zeta_acc,predict_zeta,correct_zeta,M,dt);
    ABAb_update(clt,crd,qrot,numclut,numatom);
    //    if (TERMMODE==ON) ABAb_update_Term(crd,delta_Term,numatom);

    if (TERMMODE==ON)
      ABAb_calcKineE_TermOn_new_simp(&KE,clt,crd,qvel,vel_Term,numclut,numatom);
    else
      ABAb_calcKineE(&KE,clt,crd,qvel,numclut,numatom);
    if (MODE==NVT)
      ABAbNH_calcKE_new(gzi,s,s_vel,Q_NH,KEobj,&PEv,&KEv);
    //      ABAbNH_chain_calcKE_new(zeta,zeta_vel,Q_NH,M,DOF,KBT,&PEv,&KEv);

    if (i%interval==0) {
      E_t=e_GOLM.p_MB+KE+KEv+PEv;
      p_t=e_GOLM.p_MB;
      //      E_t=e_GOLM.p_t+KE+KEv+PEv;
      //      p_t=e_GOLM.p_t;                           
      fprintf(outputfile,"%d %e %e %e %e %e %e \n",i+1,/*e_GOLM.*/p_t,KE,KEv,PEv,E_t,T);
    }

    if (i%interval2==0) {
      E_t=p_t+KE+KEv+PEv;
      fprintf(outputfile2,"E_t    = %e \n",E_t);
      fprintf(outputfile2,"KE     = %e \n",KE);
      fprintf(outputfile2,"KEv    = %e \n",KEv);
      fprintf(outputfile2,"PEv    = %e \n",PEv);
      fprintf(outputfile2,"p_tot  = %e \n",p_t);
      /************************************************************/
      /* fprintf(outputfile2,"p_nat  = %e \n",e_GOLM.p_natatt_t); */
      /* fprintf(outputfile2,"p_rep  = %e \n",e_GOLM.p_repul_t);  */
      /* fprintf(outputfile2,"p_dih  = %e \n",e_GOLM.p_d_t);	  */
      /* fprintf(outputfile2,"p_ang  = %e \n",e_GOLM.p_a_t);	  */
      /* fprintf(outputfile2,"p_bon  = %e \n",e_GOLM.p_b_t);	  */
      /************************************************************/
      fprintf(outputfile2,"T      = %e \n",T);
              
      avePE=(i*avePE+p_t)/(i+1);
      varPE=(i*varPE+p_t*p_t)/(i+1);
       
      aveKE=(i*aveKE+KE)/(i+1);
      varKE=(i*varKE+KE*KE)/(i+1);
      
      aveT=(i*aveT+T)/(i+1);
      varT=(i*varT+T*T)/(i+1);
    }
      
    if (i%interval3==0) {
      for (j=0;j<3;++j) COM[j]=0.0;
      for (j=0;j<numatom;++j) 
	for (k=0;k<3;++k) 
	  COM[k]+=mass[j]*crd[j*3+k]/summass;
      //      for (j=0;j<numatom;++j) 
      //	for (k=0;k<3;++k) 
      //	  crd[j*3+k]-=COM[k];

      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crd[j*3+k]-COM[k];
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

  logfile=efopen(logfilename,"w");
  fprintf(logfile,"initial   coordinate file = %s\n",inputfilename);
  fprintf(logfile,"reference coordinate file1 = %s\n",reffilename1);
  fprintf(logfile,"reference coordinate file2 = %s\n",reffilename2);
  fprintf(logfile,"topology             file = %s\n",parmfilename);
  fprintf(logfile,"information          file = %s\n",outputfilename);
  fprintf(logfile,"output    data       file = %s\n",outputfilename2);
  fprintf(logfile,"trajectory           file = %s\n",trjfilename);
  varPE=sqrt(varPE);
  fprintf(logfile,"Potential energy = %10.5lf kcal/mol +- %10.5lf kcal/mol\n",avePE,varPE);
  varKE=sqrt(varKE);
  fprintf(logfile,"Kinetic   energy = %10.5lf kcal/mol +- %10.5lf kcal/mol\n",aveKE,varKE);
  varT=sqrt(varT);
  fprintf(logfile,"Temperature      = %10.5lf K        +- %10.5lf K\n",aveT,varT);
  fclose(logfile);

  ABAbs_restat_write_vel_new(rstvelfilename,numclut,correct,correct_Term,correct_Term2,correct_gzi,MODE,TERMMODE);
  //  ABAbs_restat_write_vel_chain(rstvelfilename,numclut,correct,correct_Term,correct_Term2,correct_zeta,M,MODE,TERMMODE);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename refcrdfilename clustfilename parmfilename outputfilename1 outputfilenam2 trjfilename\n",progname);
}
