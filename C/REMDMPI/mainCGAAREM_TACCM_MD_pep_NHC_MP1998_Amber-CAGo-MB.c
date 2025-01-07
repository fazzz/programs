
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>


#include "REMDCGAA_TACCM_MPI_2_Amber-CAGo-MB.h"
#include "REMDCGAA_TACCM_MPI_2_Amber-CAGo-MB_2.h"
#include "REMD_functions.h"

#include "EXT_SYS.h"

#include "TACCM_CGAAMDrun_Amber_CAGo_MB_CTheta.h"

#define ON  0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d,dummy;
  int numstep=1000,interval=100;
  int numEX=1,numRE=2;
  int vMode=OFF;
  int numEXT=1;

  struct AADataforREMD_Amber_CAGo_MB AAdata;
  struct CGDataforREMD_Amber_CAGo_MB CGdata;
  struct AACGCommonDataforREMD_Amber_CAGo_MB Cdata;
  struct TACCMDataforREMD_Amber_CAGo_MB Zdata;

  /********************************************************/
  /* struct AADataforREMD_testb AAdata_test, CGdata_test; */
  /* struct AACGCommonDataforREMD_testb Cdata_test;	  */
  /* struct TACCMDataforREMD_testb Zdata_test;		  */
  /********************************************************/
  
  struct AmberParmL ap_AA,ap_CG;

  double *KZAA,*KZCG;
  double parameterCG=0.2;

  int nc=1;                          
  double T0AA,T0CG,T0Z;
  double k_B=1.98723e-3;             
  double UNITT=418.4070;             
  double KTAA,KTCG,KBTZ,tau=0.01,tau2,pi;                    

  double dt,dt2,wdt2[3],wdt4[3];

  double *refcrd1,*refcrd2,*refcrdAA1,*refcrdAA2,ep=0.3,criteria=6.5;
  double de=1.0,d_CG=1.0,d2,cutoff=6.5;
  int NCmode=3,nibnum=3,numnb,num14;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*refcrdfilename1,*refcrdfilename2,*parmfilenameAA,*parmfilenameCG;

  char *outputfilenameAAbase,*trjfilenameAAbase,outputfilenameAA[2000],trjfilenameAA[2000];
  char *outputfilenameCGbase,*trjfilenameCGbase,outputfilenameCG[2000],trjfilenameCG[2000];
  char *trjfileZbase,*trjfileThetaAAbase,*trjfileThetaCGbase,*logfilename;
  char trjfilenameZ[2000],trjfilenameThetaAA[2000],trjfilenameThetaCG[2000],logf[1000];
  FILE *inputfile,*refcrdfile1,*refcrdfile2,*parmfileAA,*parmfileCG,*logfile;

  char *acc_ratio_filename="AccRatio.txt";
  FILE *acc_ratio_file;

  char *progname;
  int opt_idx=1;

  double **acc_ratio;

  double *theta,EAA,ECG,EZ;

  int my_rank,num_procs,tag = 0;
  MPI_Status status;            

  int num_a_prot,NUM_A_PROT;
  double *ele_ele,*ALJ_s,*BLJ_s,*ele_ele_14,*ALJ_14_s,*BLJ_14_s;

  MPI_Init(&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  pi=acos(-1.0);

  struct option long_opt[] = {
    {"vMode",1,NULL,'v'},
    {"nums",1,NULL,'s'}, {"int",1,NULL,'i'},
    {"numRE",1,NULL,'n'}, {"numEX",1,NULL,'e'},
    {"numEXT",1,NULL,'E'}, {"pCG",1,NULL,'P'},
    {"tau",1,NULL,'a'}, {"dt",1,NULL,'x'},
    {"TAA",1,NULL,'T'}, {"TCG",1,NULL,'t'},  {"TZ",1,NULL,'B'},
    {"mZ",1,NULL,'m'}, {"massX",1,NULL,'X'},
    {"ep",1,NULL,'p'}, {"cutoff",1,NULL,'l'},
    {"de",1,NULL,'d'}, {"dd",1,NULL,'2'},
    {"AR",1,NULL,'A'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hv:n:E:s:e:a:i:l:x:T:t:B:m:X:p:c:E:P:d:2:A:",long_opt,&opt_idx))!=-1) {
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
    case 'd':
      d=atof(optarg);
      break;
    case '2':
      de=atof(optarg);
      break;
    case 'l':
      cutoff=atof(optarg);
      break;
    case 'x':
      dt=atof(optarg);   break;
    case 'm':
      Zdata.massZ=atof(optarg); break;
    case 'p':
      ep=atof(optarg); break;
    case 'c':
      criteria=atof(optarg); break;
    case 'E':
      numEXT=atoi(optarg);  break;
    case 'P':
      parameterCG=atof(optarg);  break;
    case 'A':
      acc_ratio_filename=optarg;  break;
    case 'h':
      USAGE(progname);  exit(1);
    default:
      USAGE(progname);   exit(1);
    }
  }

  progname=*argv;  argc-=optind;  argv+=optind;

  if (argc < 12) {
    USAGE(progname);
    exit(1);
  }
  inputfilename        = *argv;
  refcrdfilename1      = *++argv;			   
  refcrdfilename2      = *++argv;
  parmfilenameAA       = *++argv;
  outputfilenameAAbase = *++argv;
  trjfilenameAAbase    = *++argv;
  outputfilenameCGbase = *++argv;
  trjfilenameCGbase    = *++argv;
  trjfileZbase         = *++argv;
  trjfileThetaAAbase   = *++argv;
  trjfileThetaCGbase   = *++argv;
  logfilename          = *++argv;

  index_replicas=(int *)gcemalloc(sizeof(int)*numRE);
  index_parameters=(int *)gcemalloc(sizeof(int)*numRE);
  REMD_ini_purmutation_funcs(numRE);
  
  parmfileAA=efopen(parmfilenameAA,"r");
  readParmtopLb(parmfileAA,&ap_AA);
  fclose(parmfileAA);
  
  parmfileCG=efopen(parmfilenameAA,"r");
  readParmtopLb(parmfileCG,&ap_CG);
  fclose(parmfileCG);
  
  if ( ap_AA.NATOM != ap_CG.NATOM ) {
    printf("error: about num of atoms\n");
  }
  
  j=0;
  for (i=0;i<ap_CG.NATOM;++i) {
    if (strncmp(ap_CG.IGRAPH[i],"CA",2)==0) {
      ++j;
    }
  }
  Cdata.numCAatom=j;
  Cdata.numres=ap_CG.NATOM;
  Cdata.massCA=(double *)gcemalloc(sizeof(double)*Cdata.numCAatom);
  j=0;
  for (i=0;i<ap_CG.NATOM;++i) {
    if (strncmp(ap_CG.IGRAPH[i],"CA",2)==0) {
      Cdata.massCA[j]=ap_CG.AMASS[i];
      ++j;
    }
  }
  
  Cdata.numatom=ap_AA.NATOM;
  j=0;  for (i=0;i<Cdata.numatom;++i)  if (strncmp(ap_AA.IGRAPH[i],"H",1)==0)  ++j;
  Cdata.numheavyatom=Cdata.numatom-j;  Cdata.numres=ap_AA.NRES;
  
  Zdata.numZ_bon=Cdata.numCAatom-1;
  Zdata.numZ_ang=Cdata.numCAatom-2;
  Zdata.numZ_dih=Cdata.numCAatom-3;

  Zdata.pairs_dihed_AA=(int **)gcemalloc(sizeof(int *)*Zdata.numZ_dih);
  Zdata.pairs_dihed_CG=(int **)gcemalloc(sizeof(int *)*Zdata.numZ_dih);
  for (i=0;i<Zdata.numZ_dih;++i) {
     Zdata.pairs_dihed_AA[i]=(int *)gcemalloc(sizeof(int)*4);
     Zdata.pairs_dihed_CG[i]=(int *)gcemalloc(sizeof(int)*4);
  }
  Zdata.pairs_angle_AA=(int **)gcemalloc(sizeof(int *)*Zdata.numZ_ang);
  Zdata.pairs_angle_CG=(int **)gcemalloc(sizeof(int *)*Zdata.numZ_ang);
  for (i=0;i<Zdata.numZ_ang;++i) {
    Zdata.pairs_angle_AA[i]=(int *)gcemalloc(sizeof(int)*3);
    Zdata.pairs_angle_CG[i]=(int *)gcemalloc(sizeof(int)*3);
  }

  Zdata.pairs_bond_AA=(int **)gcemalloc(sizeof(int *)*Zdata.numZ_bon);
  Zdata.pairs_bond_CG=(int **)gcemalloc(sizeof(int *)*Zdata.numZ_bon);
  for (i=0;i<Zdata.numZ_bon;++i) {
    Zdata.pairs_bond_AA[i]=(int *)gcemalloc(sizeof(int)*2);
    Zdata.pairs_bond_CG[i]=(int *)gcemalloc(sizeof(int)*2);
  }

  for (i=0;i<Zdata.numZ_bon;++i) {
    Zdata.pairs_bond_CG[i][0]=i+1;
    Zdata.pairs_bond_CG[i][1]=i+2;
    Zdata.pairs_bond_AA[i][0]=PTL_canum_fromresnumb(Zdata.pairs_bond_CG[i][0],ap_AA)+1;
    Zdata.pairs_bond_AA[i][1]=PTL_canum_fromresnumb(Zdata.pairs_bond_CG[i][1],ap_AA)+1;
  }
  
  for (i=0;i<Zdata.numZ_ang;++i) {
    Zdata.pairs_angle_CG[i][0]=i+1;
    Zdata.pairs_angle_CG[i][1]=i+2;
    Zdata.pairs_angle_CG[i][2]=i+3;
    Zdata.pairs_angle_AA[i][0]=PTL_canum_fromresnumb(Zdata.pairs_angle_CG[i][0],ap_AA)+1;
    Zdata.pairs_angle_AA[i][1]=PTL_canum_fromresnumb(Zdata.pairs_angle_CG[i][1],ap_AA)+1;
    Zdata.pairs_angle_AA[i][2]=PTL_canum_fromresnumb(Zdata.pairs_angle_CG[i][2],ap_AA)+1;
  }
  
  for (i=0;i<Zdata.numZ_dih;++i) {
    Zdata.pairs_dihed_CG[i][0]=i+1;
    Zdata.pairs_dihed_CG[i][1]=i+2;
    Zdata.pairs_dihed_CG[i][2]=i+3;
    Zdata.pairs_dihed_CG[i][3]=i+4;
    Zdata.pairs_dihed_AA[i][0]=PTL_canum_fromresnumb((Zdata.pairs_dihed_CG[i][0]),ap_AA)+1;
    Zdata.pairs_dihed_AA[i][1]=PTL_canum_fromresnumb((Zdata.pairs_dihed_CG[i][1]),ap_AA)+1;
    Zdata.pairs_dihed_AA[i][2]=PTL_canum_fromresnumb((Zdata.pairs_dihed_CG[i][2]),ap_AA)+1;
    Zdata.pairs_dihed_AA[i][3]=PTL_canum_fromresnumb((Zdata.pairs_dihed_CG[i][3]),ap_AA)+1;
  }

  // check point !!

  Cdata.mass=(double *)gcemalloc(sizeof(double)*Cdata.numatom);
  for (i=0;i<Cdata.numatom;++i) Cdata.mass[i]=ap_AA.AMASS[i];
  refcrd1=(double *)gcemalloc(sizeof(double)*Cdata.numCAatom*3);
  refcrd2=(double *)gcemalloc(sizeof(double)*Cdata.numCAatom*3);
  refcrdAA1=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  refcrdAA2=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3);
  
  AAdata.crd=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3*2);
  AAdata.vel=(double *)gcemalloc(sizeof(double)*Cdata.numatom*3*2);
  CGdata.crd=(double *)gcemalloc(sizeof(double)*Cdata.numCAatom*3*2);
  CGdata.vel=(double *)gcemalloc(sizeof(double)*Cdata.numCAatom*3*2);
  
  Zdata.Z=(double *)gcemalloc(sizeof(double)*(Zdata.numZ_dih+Zdata.numZ_ang+Zdata.numZ_bon));
  Zdata.velZ=(double *)gcemalloc(sizeof(double)*(Zdata.numZ_dih+Zdata.numZ_ang+Zdata.numZ_bon));
  
  KZAA=(double *)gcemalloc(sizeof(double)*numRE);
  KZCG=(double *)gcemalloc(sizeof(double)*numRE);
  inputfile=efopen(inputfilename,"r");
  CGAAREMDreadInputs_Amber_CAGo_MB(inputfile,Cdata.numatom,numRE,my_rank,
				   AAdata.crd,AAdata.vel,CGdata.crd,CGdata.vel,
    				   KZAA,KZCG,ap_CG);
  fclose(inputfile);
  Zdata.KZAA=KZAA[my_rank];
  Zdata.KZCG=KZCG[my_rank];
  
  refcrdfile1=efopen(refcrdfilename1,"r");
  //  getline(&line,&len,refcrdfile1);
  fscanf(refcrdfile1,"%d",&dummy);
  j=0;
  for (i=0;i<Cdata.numatom;++i) {
    for (k=0;k<3;++k) fscanf(refcrdfile1,"%lf",&refcrdAA1[i*3+k]);
    if (strncmp(ap_CG.IGRAPH[i],"CA",2)==0) {
      for (k=0;k<3;++k) refcrd1[j*3+k]=refcrdAA1[i*3+k];
      ++j;
    }
  }
  fclose(refcrdfile1);

  refcrdfile2=efopen(refcrdfilename2,"r");
  //  getline(&line,&len,refcrdfile2);
  fscanf(refcrdfile2,"%d",&dummy);
  j=0;
  for (i=0;i<Cdata.numatom;++i) {
    for (k=0;k<3;++k) fscanf(refcrdfile2,"%lf",&refcrdAA2[i*3+k]);
    if (strncmp(ap_CG.IGRAPH[i],"CA",2)==0) {
      for (k=0;k<3;++k) refcrd2[j*3+k]=refcrdAA2[i*3+k];
      ++j;
    }
  }
  fclose(refcrdfile2);
  
  if ( vMode==OFF ) {
    for (i=0;i<numRE;++i) {
      MD_Generate_inivelo(AAdata.vel,Cdata.mass,Cdata.numatom,k_B*T0AA*UNITT);
      MD_Generate_inivelo(CGdata.vel,Cdata.massCA,Cdata.numCAatom,k_B*T0CG*UNITT);
      TACCM_MD_Generate_inivelo(Zdata.velZ,Zdata.massZ,
				(Zdata.numZ_dih+Zdata.numZ_ang+Zdata.numZ_bon),k_B*T0Z*UNITT);
  
      AAdata.zeta=0.0;      AAdata.V_zeta=0.0;
      CGdata.zeta=0.0;      CGdata.V_zeta=0.0;
      Zdata.zetaZ=0.0;      Zdata.V_zetaZ=0.0;
    }
  }
  TACCM_CTheta_Amber_CAGo_MB_2(AAdata.crd,Cdata.numatom,Zdata.Z,
  			       Zdata.numZ_dih,Zdata.pairs_dihed_AA,
  			       Zdata.numZ_ang,Zdata.pairs_angle_AA,
  			       Zdata.numZ_bon,Zdata.pairs_bond_AA,
  			       pi);
  
  tau=tau/2.0/pi;  tau2=tau*tau;
  KBTZ=k_B*T0Z;  Zdata.NfKTZ=((Zdata.numZ_dih+Zdata.numZ_ang+Zdata.numZ_bon)+1)*KBTZ*UNITT;
  Zdata.QZ=tau2*KBTZ*UNITT*(Zdata.numZ_dih+Zdata.numZ_ang+Zdata.numZ_bon);
  KTAA=k_B*T0AA;  AAdata.NfKT=(3.0*Cdata.numatom+1)*KTAA*UNITT;  
  AAdata.Q=tau2*KTAA*UNITT*(3.0*Cdata.numatom);
  
  KTCG=k_B*T0CG;  CGdata.NfKT=(3.0*Cdata.numCAatom+1)*KTCG*UNITT;
  CGdata.Q=tau2*KTCG*UNITT*(3.0*Cdata.numCAatom);

  ffLc_set_calcffandforce_HS(&(AAdata.e),&(AAdata.f),ap_AA,
  			     ele_ele,ALJ_s,BLJ_s,
  			     ele_ele_14,ALJ_14_s,BLJ_14_s);

  /********************************************************************************************/
  /* ele_ele=(double *)gcemalloc(sizeof(double)*(AAdata.e).parm.numnb);			      */
  /* ele_ele_14=(double *)gcemalloc(sizeof(double)*(AAdata.e).parm.num14);		      */
  /* 											      */
  /* ALJ_s=(double *)gcemalloc(sizeof(double)*(AAdata.e).parm.numnb);			      */
  /* BLJ_s=(double *)gcemalloc(sizeof(double)*(AAdata.e).parm.numnb);			      */
  /* ALJ_14_s=(double *)gcemalloc(sizeof(double)*(AAdata.e).parm.num14);		      */
  /* BLJ_14_s=(double *)gcemalloc(sizeof(double)*(AAdata.e).parm.num14);		      */
  /* 											      */
  /* for (i=0;i<(AAdata.e).parm.numnb;++i) {						      */
  /*   num_a_prot=(AAdata.e).parm.indexnb[i*2];						      */
  /*   NUM_A_PROT=(AAdata.e).parm.indexnb[i*2+1];					      */
  /* 											      */
  /*   ele_ele[i]=(AAdata.e).parm.ele[num_a_prot]*(AAdata.e).parm.ele[NUM_A_PROT];	      */
  /*   ALJ_s[i]=(AAdata.e).parm.ALJ[num_a_prot*num_a_prot+NUM_A_PROT];			      */
  /*   BLJ_s[i]=(AAdata.e).parm.BLJ[num_a_prot*num_a_prot+NUM_A_PROT];			      */
  /* }											      */
  /* 											      */
  /* for (i=0;i<(AAdata.e).parm.num14;++i) {						      */
  /*   num_a_prot=(AAdata.e).parm.index14[i*2];						      */
  /*   NUM_A_PROT=(AAdata.e).parm.index14[i*2+1];					      */
  /* 											      */
  /*   ele_ele_14[i]=1.0/1.2*(AAdata.e).parm.ele[num_a_prot]*(AAdata.e).parm.ele[NUM_A_PROT]; */
  /*   ALJ_14_s[i]=0.5*(AAdata.e).parm.ALJ[num_a_prot*num_a_prot+NUM_A_PROT];		      */
  /*   BLJ_14_s[i]=0.5*(AAdata.e).parm.BLJ[num_a_prot*num_a_prot+NUM_A_PROT];		      */
  /* }											      */
  /********************************************************************************************/

  GOLM_Clementi_MB_ff_set_calcff_wmCutOff(&(CGdata.e_GOLM),refcrd1,refcrd2,refcrdAA1,refcrdAA2,
					  Cdata.numCAatom,Cdata.numatom,ep,cutoff);
  
  d2=d*d;
  d2=d2*k_B;
  de=de*k_B;
  
  MDb_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);
  
  sprintf(outputfilenameAA,"%s_%d",outputfilenameAAbase,my_rank+1);
  sprintf(trjfilenameAA,"%s_%d",trjfilenameAAbase,my_rank+1);
  
  sprintf(outputfilenameCG,"%s_%d",outputfilenameCGbase,my_rank+1);
  sprintf(trjfilenameCG,"%s_%d",trjfilenameCGbase,my_rank+1);
  
  sprintf(trjfilenameZ,"%s_%d",trjfileZbase,my_rank+1);
  sprintf(trjfilenameThetaAA,"%s_%d",trjfileThetaAAbase,my_rank+1);
  sprintf(trjfilenameThetaCG,"%s_%d",trjfileThetaCGbase,my_rank+1);
  
  sprintf(logf,"%s_%d_ex.log",logfilename,my_rank+1);
  
  myncL_create_def_MCD(trjfilenameAA,Cdata.numatom,&(AAdata.nc_id_MCD));
  AAdata.outputfile=efopen(outputfilenameAA,"w");
  
  myncL_create_def_MCD(trjfilenameCG,Cdata.numatom,&(CGdata.nc_id_MCD));
  CGdata.outputfile=efopen(outputfilenameCG,"w");
  
  Zdata.trjfileZ=efopen(trjfilenameZ,"w");
  Zdata.trjfilThetaAA=efopen(trjfilenameThetaAA,"w");
  Zdata.trjfilThetaCG=efopen(trjfilenameThetaCG,"w");

  MPI_Barrier(MPI_COMM_WORLD);
  
  if ( num_procs != numRE ) {    printf("condition error\n");    exit(1);  }
  
  logfile=efopen(logf,"w");
  acc_ratio=MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_Amber_CAGo_MB_2(my_rank,num_procs,tag,&status,
								 numRE,numEX,KZAA,KZCG,
								 AAdata,
								 CGdata,
								 Zdata,
								 Cdata,
								 ap_AA,ap_CG,
								 de,d2,
								 T0AA,T0CG,T0Z,
								 numstep,interval,
								 dt,dt2,wdt2,wdt4,nc,
								 UNITT,k_B,tau,pi,
								 parameterCG,logfile,
								 ele_ele,ALJ_s,BLJ_s,
								 ele_ele_14,ALJ_14_s,BLJ_14_s);
  
  fclose(logfile);
  fclose(AAdata.outputfile);
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
  printf("%s [-h] inputfilename refcrdfilename1 refcrdfilename2 TACCMfilename parmfilename outputfilenameAAbase trjfilenameAAbase outputfilenameCGbase trjfilenameCGbase trjfileZbase trjfileThetaZ\n",progname);
}
