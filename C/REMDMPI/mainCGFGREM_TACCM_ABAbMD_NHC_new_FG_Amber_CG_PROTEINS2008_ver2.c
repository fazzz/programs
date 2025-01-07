
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>
#include "mpi.h"

#include "REMDCGAA_TACCM_MPI_2_Amber_PROTEINS2008.h"
#include "REMDCGAA_TACCM_MPI_2_Amber_PROTEINS2008_ver2.h"
#include "REMDCGAA_TACCM_MPI_2_Amber_PROTEINS2008_Amber_hybrid.h"

#include "TACCM_CGFGABAbMDrun_Amber_PROTEINS2008.h"
#include "TACCM_CGFGABAbMDrun_Amber_PROTEINS2008_ver2.h"

#include "REMD_functions.h"

#include "ABAb.h"

#define ON  0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d;
  int numstep=1000,interval=100;
  int numEX=1,numRE=2;
  int vMode=OFF;

  int equflag=OFF;
  int numstepequ=0,intervalequ=1;

  struct AADataforREMD_Amber FGdata;
  struct CGDataforREMD_PROTEINS2008 CGdata;
  struct AACGCommonDataforREMD_A_P2008 Cdata;
  struct TACCMDataforREMD_A_P2008_ver2 Zdata;
  double *KZAA,*KZCG;
  double *massZs;

  struct potential e;
  struct force f;

  int nc=1;                          
  double T0AA,T0CG,T0Z;
  double k_B=1.98723e-3;             
  double UNITT=418.4070;             
  double KTAA,KTCG,KTZ,tau=0.01,tau2,pi;                    

  double dt,dt2;

  double *refcrd,ep=0.3,criteria=6.5;
  int NCmode=3,nibnum=3,numnb,num14;
  double EAAm_equ,ECGm_equ,EZm_equ;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*refcrdfilename,*parmfilename,*TACCMfilename;

  char *outputfilenameAAbase,*trjfilenameAAbase,outputfilenameAA[2000],trjfilenameAA[2000];
  char *outputfilenameCGbase,*trjfilenameCGbase,outputfilenameCG[2000],trjfilenameCG[2000];
  char *trjfilenameZcrdbase,trjfilenameZcrd[2000];
  char *trjfileZbase,*trjfileThetaAAbase,*trjfileThetaCGbase,*logfilename;
  char trjfilenameZ[2000],trjfilenameThetaAA[2000],trjfilenameThetaCG[2000],logf[1000];
  char *clustfilename;
  FILE *inputfile,*refcrdfile,*parmfile,*TACCMfile,*logfile;
  FILE *clustfile;

  char *acc_ratio_filename="AccRatio.txt";
  FILE *acc_ratio_file;

  char *progname;
  int opt_idx=1;

  double **acc_ratio;

  int my_rank,num_procs,tag = 0;
  MPI_Status status;            

  MPI_Init(&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  pi=acos(-1.0);

  struct option long_opt[] = {
    {"vMode",1,NULL,'v'},
    {"nums",1,NULL,'s'}, {"int",1,NULL,'i'},
    {"numRE",1,NULL,'n'}, {"numEX",1,NULL,'e'},
    {"tau",1,NULL,'a'}, {"dt",1,NULL,'x'},
    {"TAA",1,NULL,'T'}, {"TCG",1,NULL,'t'},  {"TZ",1,NULL,'B'},
    {"mZ",1,NULL,'m'}, {"massX",1,NULL,'X'},
    {"ep",1,NULL,'p'}, {"cutoff",1,NULL,'c'},
    {"AR",1,NULL,'A'},
    {"equ",1,NULL,'E'}, {"intequ",1,NULL,'I'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hv:n:s:e:a:i:x:T:t:B:m:X:p:c:A:E:I:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 's':
      numstep=atoi(optarg); break;
    case 'i':
      interval=atoi(optarg);  break;
    case 'e':
      numEX=atoi(optarg); break;
    case 'n':
      numRE=atoi(optarg); break;
    case 'T':
      T0AA=atof(optarg);  break;
    case 't':
      T0CG=atof(optarg);  break;
    case 'B':
      T0Z=atof(optarg);break;
    case 'v':
      vMode=ON;        break;
    case 'a':
      tau=atof(optarg); break;
    case 'x':
      dt=atof(optarg);   break;
    case 'm':
      Zdata.massZ=atof(optarg);  break;
    case 'p':
      ep=atof(optarg); break;
    case 'c':
      criteria=atof(optarg); break;
    case 'A':
      acc_ratio_filename=optarg;  break;
    case 'E':
      equflag=ON;
      numstepequ=atoi(optarg);  break;
    case 'I':
      intervalequ=atoi(optarg);  break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);   exit(1);
    }
  }

  progname=*argv;  argc-=optind;  argv+=optind;

  if (argc < 14) {
    USAGE(progname);
    exit(1);
  }
  inputfilename        = *argv;
  refcrdfilename       = *++argv;			   
  clustfilename        = *++argv;			   
  TACCMfilename        = *++argv;			 
  parmfilename         = *++argv;
  outputfilenameAAbase = *++argv;
  trjfilenameAAbase    = *++argv;
  outputfilenameCGbase = *++argv;
  trjfilenameCGbase    = *++argv;
  trjfilenameZcrdbase  = *++argv;
  trjfileZbase         = *++argv;
  trjfileThetaAAbase   = *++argv;
  trjfileThetaCGbase   = *++argv;
  logfilename          = *++argv;

  //  printf("%d:165\n",my_rank);

  dt2=dt*dt;

  index_replicas=(int *)gcemalloc(sizeof(int)*numRE);
  index_parameters=(int *)gcemalloc(sizeof(int)*numRE);
  REMD_ini_purmutation_funcs(numRE);

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 

  (Cdata).numatom=AP.NATOM;
  j=0;  for (i=0;i<Cdata.numatom;++i)  if (strncmp(AP.IGRAPH[i],"H",1)==0)  ++j;
  Cdata.numheavyatom=Cdata.numatom-j;  Cdata.numres=AP.NRES;

  Cdata.mass=(double *)gcemalloc(sizeof(double)*Cdata.numatom);
  for (i=0;i<Cdata.numatom;++i) Cdata.mass[i]=AP.AMASS[i];
  refcrd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);

  TACCMfile=efopen(TACCMfilename,"r");
  fscanf(TACCMfile,"%d",&Zdata.numZ);
  Zdata.pairs=(int **)gcemalloc(sizeof(int *)*Zdata.numZ);
  for (i=0;i<Zdata.numZ;++i) Zdata.pairs[i]=(int *)gcemalloc(sizeof(int)*5);
  for (i=0;i<Zdata.numZ;++i) {
    for (j=0;j<4;++j) fscanf(TACCMfile,"%d",&(Zdata.pairs[i][j]));
    fscanf(TACCMfile,"%d",&(Zdata.pairs[i][j]));
  }
  fclose(TACCMfile);

  refcrdfile=efopen(refcrdfilename,"r");
  getline(&line,&len,refcrdfile);
  fscanf(refcrdfile,"%d",&d);
  for (i=0;i<Cdata.numatom;++i) for (j=0;j<3;++j) fscanf(refcrdfile,"%lf",&refcrd[i*3+j]);
  fclose(refcrdfile);

  clustfile=efopen(clustfilename,"r");
  (FGdata.clt)=ABAbp_clustscan(clustfile,&(Cdata.numclut));
  CGdata.clt=(CLTb *)gcemalloc(sizeof(CLTb)*(Cdata.numclut));
  Zdata.clt=(CLTb *)gcemalloc(sizeof(CLTb)*(Cdata.numclut));

  for(i=0;i<Cdata.numclut;++i) {
    CGdata.clt[i].origin_atom_a=FGdata.clt[i].origin_atom_a;
    Zdata.clt[i].origin_atom_a=FGdata.clt[i].origin_atom_a;
  }
  for(i=0;i<Cdata.numclut;++i) {
    CGdata.clt[i].terminal=FGdata.clt[i].terminal;
    Zdata.clt[i].terminal=FGdata.clt[i].terminal;
  }
  for(i=0;i<Cdata.numclut;++i) {
    CGdata.clt[i].num_atom_clust=FGdata.clt[i].num_atom_clust;
    Zdata.clt[i].num_atom_clust=FGdata.clt[i].num_atom_clust;
  }
  for(i=0;i<Cdata.numclut;++i) {
    CGdata.clt[i].num_branch=FGdata.clt[i].num_branch;
    CGdata.clt[i].terminal_atom_a=(int *)gcemalloc(sizeof(int)*CGdata.clt[i].num_branch);
    CGdata.clt[i].nNumClutOfChild=(int *)gcemalloc(sizeof(int)*CGdata.clt[i].num_branch);

    Zdata.clt[i].num_branch=FGdata.clt[i].num_branch;
    Zdata.clt[i].terminal_atom_a=(int *)gcemalloc(sizeof(int)*CGdata.clt[i].num_branch);
    Zdata.clt[i].nNumClutOfChild=(int *)gcemalloc(sizeof(int)*CGdata.clt[i].num_branch);
  }
  for(i=0;i<Cdata.numclut;++i) {
    for(j=0;j<CGdata.clt[i].num_branch;++j) {
      CGdata.clt[i].terminal_atom_a[j]=FGdata.clt[i].terminal_atom_a[j];
      Zdata.clt[i].terminal_atom_a[j]=FGdata.clt[i].terminal_atom_a[j];
    }
  }
  for (i=0;i<Cdata.numclut;++i) {
    CGdata.clt[i].nNumClutOfParent=FGdata.clt[i].nNumClutOfParent;
    Zdata.clt[i].nNumClutOfParent=FGdata.clt[i].nNumClutOfParent;
  }
  for (i=0;i<Cdata.numclut;++i) {
    for(j=0;j<CGdata.clt[i].num_branch;++j) {
      CGdata.clt[i].nNumClutOfChild[j]=FGdata.clt[i].nNumClutOfChild[j];
      Zdata.clt[i].nNumClutOfChild[j]=FGdata.clt[i].nNumClutOfChild[j];
    }
  }

  for(i=0;i<Cdata.numclut;++i) {
    CGdata.clt[i].frc_e=(double *)gcemalloc(sizeof(double)*CGdata.clt[i].num_atom_clust*3);
    CGdata.clt[i].frc_LJ=(double *)gcemalloc(sizeof(double)*CGdata.clt[i].num_atom_clust*3);
    CGdata.clt[i].frc_1_4_e=(double *)gcemalloc(sizeof(double)*CGdata.clt[i].num_atom_clust*3);
    CGdata.clt[i].frc_1_4_LJ=(double *)gcemalloc(sizeof(double)*CGdata.clt[i].num_atom_clust*3);

    Zdata.clt[i].frc_e=(double *)gcemalloc(sizeof(double)*CGdata.clt[i].num_atom_clust*3);
    Zdata.clt[i].frc_LJ=(double *)gcemalloc(sizeof(double)*CGdata.clt[i].num_atom_clust*3);
    Zdata.clt[i].frc_1_4_e=(double *)gcemalloc(sizeof(double)*CGdata.clt[i].num_atom_clust*3);
    Zdata.clt[i].frc_1_4_LJ=(double *)gcemalloc(sizeof(double)*CGdata.clt[i].num_atom_clust*3);
  }
  fclose(clustfile);

  FGdata.crd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  FGdata.q=(double *)gcemalloc(sizeof(double)*Cdata.numclut);
  FGdata.qvel=(double *)gcemalloc(sizeof(double)*Cdata.numclut);
  FGdata.predict=(double *)gcemalloc(sizeof(double)*Cdata.numclut*6);
  FGdata.correct=(double *)gcemalloc(sizeof(double)*Cdata.numclut*6);

  FGdata.vel_Term=(double *)gcemalloc(sizeof(double)*6);
  FGdata.predict_Term=(double **)gcemalloc(sizeof(double *)*6);
  FGdata.predict_Term2=(double **)gcemalloc(sizeof(double *)*6);
  FGdata.correct_Term=(double **)gcemalloc(sizeof(double *)*6);
  FGdata.correct_Term2=(double **)gcemalloc(sizeof(double *)*6);
  for (i=0;i<6;++i) {
    FGdata.predict_Term[i]=(double *)gcemalloc(sizeof(double)*6);
    FGdata.predict_Term2[i]=(double *)gcemalloc(sizeof(double)*6);
    FGdata.correct_Term[i]=(double *)gcemalloc(sizeof(double)*6);
    FGdata.correct_Term2[i]=(double *)gcemalloc(sizeof(double)*6);
  }

  for (i=0;i<6;++i) FGdata.vel_Term[i]=0.0;
  for (i=0;i<(Cdata).numclut;++i) FGdata.qvel[i]=0.00;

  CGdata.crd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  CGdata.q=(double *)gcemalloc(sizeof(double)*Cdata.numclut);
  CGdata.qvel=(double *)gcemalloc(sizeof(double)*Cdata.numclut);
  CGdata.predict=(double *)gcemalloc(sizeof(double)*Cdata.numclut*6);
  CGdata.correct=(double *)gcemalloc(sizeof(double)*Cdata.numclut*6);

  CGdata.vel_Term=(double *)gcemalloc(sizeof(double)*6);
  CGdata.predict_Term=(double **)gcemalloc(sizeof(double *)*6);
  CGdata.predict_Term2=(double **)gcemalloc(sizeof(double *)*6);
  CGdata.correct_Term=(double **)gcemalloc(sizeof(double *)*6);
  CGdata.correct_Term2=(double **)gcemalloc(sizeof(double *)*6);
  for (i=0;i<6;++i) {
    CGdata.predict_Term[i]=(double *)gcemalloc(sizeof(double)*6);
    CGdata.predict_Term2[i]=(double *)gcemalloc(sizeof(double)*6);
    CGdata.correct_Term[i]=(double *)gcemalloc(sizeof(double)*6);
    CGdata.correct_Term2[i]=(double *)gcemalloc(sizeof(double)*6);
  }

  for (i=0;i<6;++i) CGdata.vel_Term[i]=0.0;
  for (i=0;i<(Cdata).numclut;++i) CGdata.qvel[i]=0.00;

  Zdata.crd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  Zdata.Z=(double *)gcemalloc(sizeof(double)*Zdata.numZ);
  Zdata.q=(double *)gcemalloc(sizeof(double)*Cdata.numclut);
  Zdata.qvel=(double *)gcemalloc(sizeof(double)*Cdata.numclut);
  Zdata.predict=(double *)gcemalloc(sizeof(double)*Cdata.numclut*6);
  Zdata.correct=(double *)gcemalloc(sizeof(double)*Cdata.numclut*6);

  Zdata.vel_Term=(double *)gcemalloc(sizeof(double)*6);
  Zdata.predict_Term=(double **)gcemalloc(sizeof(double *)*6);
  Zdata.predict_Term2=(double **)gcemalloc(sizeof(double *)*6);
  Zdata.correct_Term=(double **)gcemalloc(sizeof(double *)*6);
  Zdata.correct_Term2=(double **)gcemalloc(sizeof(double *)*6);
  for (i=0;i<6;++i) {
    Zdata.predict_Term[i]=(double *)gcemalloc(sizeof(double)*6);
    Zdata.predict_Term2[i]=(double *)gcemalloc(sizeof(double)*6);
    Zdata.correct_Term[i]=(double *)gcemalloc(sizeof(double)*6);
    Zdata.correct_Term2[i]=(double *)gcemalloc(sizeof(double)*6);
  }

  for (i=0;i<6;++i) Zdata.vel_Term[i]=0.0;
  for (i=0;i<(Cdata).numclut;++i) Zdata.qvel[i]=0.00;

  KZAA=(double *)gcemalloc(sizeof(double)*numRE);
  KZCG=(double *)gcemalloc(sizeof(double)*numRE);
  inputfile=efopen(inputfilename,"r");
  CGAAREMDreadInputs_Amber_PROTEINS2008_Amber_hybrid(inputfile,Cdata.numatom,numRE,my_rank,
						     FGdata.crd,FGdata.vel,CGdata.crd,CGdata.vel,
						     KZAA,KZCG);
  fclose(inputfile);

  for (i=0;i<Cdata.numatom;++i)
    for (j=0;j<3;++j) 
      Zdata.crd[i*3+j]=FGdata.crd[i*3+j];

  Zdata.KZAA=KZAA[my_rank];
  Zdata.KZCG=KZCG[my_rank];

  refcrdfile=efopen(refcrdfilename,"r");
  getline(&line,&len,refcrdfile);  fscanf(refcrdfile,"%d",&d);
  for (i=0;i<Cdata.numatom;++i) for (j=0;j<3;++j) fscanf(refcrdfile,"%lf",&refcrd[i*3+j]);
  fclose(refcrdfile);

  //  TACCM_CTheta(FGdata.crd,Cdata.numatom,Zdata.ZTheta,Zdata.numZ,Zdata.pairs,pi);
  /*************************************************/
  /* for (i=0;i<Zdata.numZ;++i) {		   */
  /*   Zdata.correct[i][0]=Zdata.Z[i];		   */
  /*   //    Zdata.correct[i][1]=dt*Zdata.velZ[i]; */
  /* }						   */
  /*************************************************/

  Cdata.numclutparent=(int *)gcemalloc(sizeof(int)*(Cdata).numclut);
  Cdata.terminal=(int *)gcemalloc(sizeof(int)*(Cdata).numclut);
  Cdata.origin=(int *)gcemalloc(sizeof(int)*(Cdata).numclut);
  for (i=0;i<(Cdata).numclut;++i) {
    Cdata.numclutparent[i]=(FGdata).clt[i].nNumClutOfParent;
    Cdata.terminal[i]=(FGdata).clt[i].terminal_atom_a[0];
    Cdata.origin[i]=(FGdata).clt[i].origin_atom_a;

    Cdata.numclutparent[i]=(CGdata).clt[i].nNumClutOfParent;
    Cdata.terminal[i]=(CGdata).clt[i].terminal_atom_a[0];
    Cdata.origin[i]=(CGdata).clt[i].origin_atom_a;
  }

  (FGdata).clt[0].join=0;
  for(i=1; i<(Cdata).numclut; ++i) {
    ABAb_setJoin((FGdata).clt,i);
  }
  ABAbs_local_reference((FGdata).clt,(Cdata).numclut,(Cdata).numatom,(FGdata).crd);
  ABAbs_trans_Matrix((FGdata).clt,(Cdata).numclut,(Cdata).numatom,(FGdata).crd);
  ABAbs_inertia_matrix((FGdata).clt,(Cdata).numclut,(Cdata).numatom,(FGdata).crd,(Cdata).mass);

  (CGdata).clt[0].join=0;
  for(i=1; i<(Cdata).numclut; ++i) {
    ABAb_setJoin((CGdata).clt,i);
  }
  ABAbs_local_reference((CGdata).clt,(Cdata).numclut,(Cdata).numatom,(CGdata).crd);
  ABAbs_trans_Matrix((CGdata).clt,(Cdata).numclut,(Cdata).numatom,(CGdata).crd);
  ABAbs_inertia_matrix((CGdata).clt,(Cdata).numclut,(Cdata).numatom,(CGdata).crd,(Cdata).mass);

  massZs=(double *)gcemalloc(sizeof(double)*(Cdata).numatom);
  for(i=1; i<(Cdata).numatom; ++i) {
    massZs[i]=(Zdata).massZ;
  }  
  (Zdata).clt[0].join=0;
  for(i=1; i<(Cdata).numclut; ++i) {
    ABAb_setJoin((Zdata).clt,i);
  }
  ABAbs_local_reference((Zdata).clt,(Cdata).numclut,(Cdata).numatom,(Zdata).crd);
  ABAbs_trans_Matrix((Zdata).clt,(Cdata).numclut,(Cdata).numatom,(Zdata).crd);
  ABAbs_inertia_matrix((Zdata).clt,(Cdata).numclut,(Cdata).numatom,(Zdata).crd,massZs);

  if ( vMode==OFF ) {
    //    MD_Generate_inivelo(FGdata.vel,Cdata.mass,Cdata.numatom,k_B*T0AA*UNITT);
    //    MD_Generate_inivelo(CGdata.vel,Cdata.mass,Cdata.numatom,k_B*T0CG*UNITT);
    /*******************************************************************************/
    /* TACCM_MD_Generate_inivelo(Zdata.velZ,Zdata.massZ,Zdata.numZ,k_B*T0Z*UNITT); */
    /* for (i=0;i<Zdata.numZ;++i) {						   */
    /*   Zdata.correct[i][1]=dt*Zdata.velZ[i];					   */
    /* }									   */
    /*******************************************************************************/

    //    FGdata.zeta=0.0;      FGdata.V_zeta=0.0;
    //    CGdata.zeta=0.0;      CGdata.V_zeta=0.0;
    //    Zdata.zetaZ=0.0;      Zdata.V_zetaZ=0.0;
  }  
  //  TACCM_CTheta(FGdata.crd,Cdata.numatom,Zdata.Z,Zdata.numZ,Zdata.pairs,pi);

  (Cdata).DOF=((Cdata).numclut-1)+6;
  //  tau=tau/2.0/pi;  tau2=tau*tau;  
  //  (Zdata).KEo=k_B*T0Z;  Zdata.NfKTZ=(Zdata.numZ+1)*(Zdata).KEo*UNITT;
  //  Zdata.QZ=tau2*(Zdata).KEo*UNITT*Zdata.numZ;
  //  (FGdata).KEo=k_B*T0AA;  FGdata.NfKT=((Cdata).DOF+1)*(FGdata).KEo*UNITT;  
  //  FGdata.Q=tau2*(FGdata).KEo*UNITT*((Cdata).DOF);
  //  (CGdata).KEo=k_B*T0CG;  CGdata.NfKT=((Cdata).DOF+1)*(CGdata).KEo*UNITT;  
  //  CGdata.Q=tau2*(CGdata).KEo*UNITT*((Cdata).DOF);

  (FGdata).KEo=0.5*(Cdata).DOF*k_B*T0AA;
  (CGdata).KEo=0.5*(Cdata).DOF*k_B*T0CG;
  (Zdata).KEo=0.5*(Cdata).DOF*k_B*T0Z;
  //  (Zdata).KEo=0.5*(Zdata).numZ*k_B*T0Z;

  ffL_set_calcffandforce(&(FGdata.e),&(FGdata.f));

  ffL_set_calcffandforce(&e,&f);
  ffL_set_non_bonding_index_1(&numnb,&num14);
  e.parm.numnb=numnb;
  e.parm.num14=num14;
  e.parm.indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  e.parm.index14=(int *)gcemalloc(sizeof(int)*num14*2);
  ffL_set_non_bonding_index_2(e.parm.indexnb,e.parm.index14);

  GOLMAA_PROTEINS2008_ff_set_calcff_b(&(CGdata.e),refcrd,Cdata.numatom,Cdata.numres,
  				      /*FGdata.*/e.parm.indexnb,/*FGdata.*/e.parm.numnb,ep,nibnum,criteria);

  ABAbNH_set_new((FGdata).s,(FGdata).s_vel,(FGdata).gzi,
		 (FGdata).predict_gzi,(FGdata).correct_gzi,
		 (FGdata).predict_s,(FGdata).correct_s,tau,&tau2,
		 &((FGdata).Q),KTAA,dt);
  ABAb_integ_set((FGdata).q,(FGdata).qvel,
		 (FGdata).predict,(FGdata).correct,
		 (Cdata).numclut,dt);

  ABAbNH_set_new((CGdata).s,(CGdata).s_vel,(CGdata).gzi,
		 (CGdata).predict_gzi,(CGdata).correct_gzi,
		 (CGdata).predict_s,(CGdata).correct_s,tau,&tau2,
		 &((CGdata).Q),KTCG,dt);
  ABAb_integ_set((CGdata).q,(CGdata).qvel,
		 (CGdata).predict,(CGdata).correct,
		 (Cdata).numclut,dt);

  ABAbNH_set_new((Zdata).s,(Zdata).s_vel,(Zdata).gzi,
		 (Zdata).predict_gzi,(Zdata).correct_gzi,
		 (Zdata).predict_s,(Zdata).correct_s,tau,&tau2,
		 &((Zdata).Q),KTZ,dt);
  ABAb_integ_set((Zdata).q,(Zdata).qvel,
		 (Zdata).predict,(Zdata).correct,
		 (Cdata).numclut,dt);

  //  printf("%d:380\n",my_rank);

  sprintf(outputfilenameAA,"%s_%d",outputfilenameAAbase,my_rank+1);
  sprintf(trjfilenameAA,"%s_%d",trjfilenameAAbase,my_rank+1);

  sprintf(outputfilenameCG,"%s_%d",outputfilenameCGbase,my_rank+1);
  sprintf(trjfilenameCG,"%s_%d",trjfilenameCGbase,my_rank+1);

  sprintf(trjfilenameZcrd,"%s_%d",trjfilenameZcrdbase,my_rank+1);

  sprintf(trjfilenameZ,"%s_%d",trjfileZbase,my_rank+1);
  sprintf(trjfilenameThetaAA,"%s_%d",trjfileThetaAAbase,my_rank+1);
  sprintf(trjfilenameThetaCG,"%s_%d",trjfileThetaCGbase,my_rank+1);

  sprintf(logf,"%s_%d_ex.log",logfilename,my_rank+1);

  myncL_create_def_MCD(trjfilenameAA,Cdata.numatom,&(FGdata.nc_id_MCD));
  FGdata.outputfile=efopen(outputfilenameAA,"w");

  myncL_create_def_MCD(trjfilenameCG,Cdata.numatom,&(CGdata.nc_id_MCD));
  CGdata.outputfile=efopen(outputfilenameCG,"w");

  myncL_create_def_MCD(trjfilenameZcrd,Cdata.numatom,&(Zdata.nc_id_MCD));

  Zdata.trjfileZ=efopen(trjfilenameZ,"w");
  Zdata.trjfilThetaAA=efopen(trjfilenameThetaAA,"w");
  Zdata.trjfilThetaCG=efopen(trjfilenameThetaCG,"w");

  MPI_Barrier(MPI_COMM_WORLD);

  if ( num_procs != numRE ) {    printf("condition error\n");    exit(1);  }

  //  printf("406\n");
  logfile=efopen(logf,"w");

  if (equflag==ON) {
    //    printf("409\n");

    runTACCM_CGFG_ABAbMD_NH_new_Amber_PROTEINS2008_ver2(// AA ////////////////////////////////////////////////////
							FGdata.crd,FGdata.q,FGdata.qvel,
							FGdata.predict,FGdata.correct,
							FGdata.s,FGdata.s_vel,FGdata.predict_s,FGdata.correct_s,
							FGdata.gzi,FGdata.gzi_vel,
							FGdata.predict_gzi,FGdata.correct_gzi,
							Cdata.numclutparent,Cdata.terminal,Cdata.origin,
							FGdata.vel_Term,FGdata.predict_Term,FGdata.predict_Term2,
							FGdata.correct_Term,FGdata.correct_Term2,
							FGdata.clt,FGdata.Q,
							FGdata.e,FGdata.f,FGdata.T,T0AA,FGdata.KEo,
							FGdata.avePE,FGdata.aveKE,FGdata.aveT,
							FGdata.varPE,FGdata.varKE,FGdata.varT,
							FGdata.nc_id_MCD,FGdata.outputfile,
							// CG ////////////////////////////////////////////////////
							CGdata.crd,CGdata.q,CGdata.qvel,
							CGdata.predict,CGdata.correct,
							CGdata.s,CGdata.s_vel,CGdata.predict_s,CGdata.correct_s,
							CGdata.gzi,CGdata.gzi_vel,
							CGdata.predict_gzi,CGdata.correct_gzi,
							Cdata.numclutparent,Cdata.terminal,Cdata.origin,
							CGdata.vel_Term,CGdata.predict_Term,CGdata.predict_Term2,
							CGdata.correct_Term,CGdata.correct_Term2,
							CGdata.clt,CGdata.Q,
							CGdata.e,CGdata.T,T0CG,CGdata.KEo,
							CGdata.avePE,CGdata.aveKE,CGdata.aveT,
							CGdata.varPE,CGdata.varKE,CGdata.varT,
							CGdata.nc_id_MCD,CGdata.outputfile,
							// Z  /////////////////////////////////////////////////
							Zdata.crd,Zdata.q,Zdata.qvel,
							Zdata.predict,Zdata.correct,
							Zdata.s,Zdata.s_vel,Zdata.predict_s,Zdata.correct_s,
							Zdata.gzi,Zdata.gzi_vel,
							Zdata.predict_gzi,Zdata.correct_gzi,
							Cdata.numclutparent,Cdata.terminal,Cdata.origin,
							Zdata.vel_Term,Zdata.predict_Term,Zdata.predict_Term2,
							Zdata.correct_Term,Zdata.correct_Term2,
							Zdata.clt,Zdata.Q,
							Zdata.T,T0Z,Zdata.KEo,
							Zdata.Z,Zdata.numZ,Zdata.massZ,
							Zdata.KZAA,Zdata.KZCG,Zdata.pairs,
							Zdata.avePE,Zdata.aveKE,Zdata.aveT,
							Zdata.varPE,Zdata.varKE,Zdata.varT, 
							Zdata.nc_id_MCD,
							Zdata.trjfileZ,Zdata.trjfilThetaAA,Zdata.trjfilThetaCG,
							// CM  ///////////////////////////////////////////////////
							Cdata.mass,(Cdata.numatom),(Cdata.numclut),(Cdata.DOF),
							numstepequ,intervalequ,&l,dt,tau,tau2,UNITT,k_B,pi,
							&EAAm_equ,&ECGm_equ,&EZm_equ);

    //    printf("451\n");

    //    FGdata.zeta=0.0; FGdata.V_zeta=0.0;
    //    CGdata.zeta=0.0; CGdata.V_zeta=0.0;
    //    Zdata.zetaZ=0.0; Zdata.V_zetaZ=0.0;

  }

  acc_ratio=MPI_CGFGTREM_TACCM_ABAbMD_NH_new_Amber_PROTEINS2008_ver2
    (my_rank,num_procs,tag,&status,
     numRE,numEX,KZAA,KZCG,
     FGdata,CGdata,Zdata,Cdata,
     T0AA,T0CG,T0Z,
     numstep,interval,dt,tau,tau2,
     UNITT,k_B,pi,logfile);
  
  fclose(logfile);
  fclose(FGdata.outputfile);  
  fclose(CGdata.outputfile);  
  fclose(Zdata.trjfileZ);  
  fclose(Zdata.trjfilThetaAA); 
  fclose(Zdata.trjfilThetaCG);

  if ( my_rank==0 ) {
    acc_ratio_file=efopen(acc_ratio_filename,"w");
    for (i=0;i<numRE;++i) {
      if (i!=numRE-1) fprintf(acc_ratio_file,"%d-%d %8.4lf\n",i,i+1,acc_ratio[i][i+1]);
      else fprintf(acc_ratio_file,"%d-0 %8.4lf\n",i,acc_ratio[i][0]);
    }
    fclose(acc_ratio_file);
  }
  
  MPI_Finalize();

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename refcrdfilename TACCMfilename parmfilename outputfilenameAAbase trjfilenameAAbase outputfilenameCGbase trjfilenameCGbase trjfileZbase trjfileThetaZ\n",progname);
}
