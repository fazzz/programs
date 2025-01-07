
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "MS_TACCM_REMD.h"

#include "REMD_functions.h"

#define ON  0
#define OFF 1

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numstep=1000,interval=100;
  int numEX=1,numRE=2;
  int vMode=OFF;

  int equflag=OFF;
  int numstepequ=0,intervalequ=1;

  int numAA=1,numCG=1;

  struct AADataMSREMD_Amber        *FGdata;
  struct CGDataMSREMD_PROTEINS2008 *CGdata;

  struct ZDataMSREMD_A_P2008_pca Zdata;

  double **KZAA,**KZCG;

  struct AmberParmL *ap,*ap_CG;
  struct potential *e;
  struct force *f;

  int nc=1;                          
  double *T0AA,*T0CG,T0Z;
  double k_B=1.98723e-3;             
  double UNITT=418.4070;             
  double *KTAA,*KTCG,KBTZ,tau=0.01,tau2,pi;                    

  double dt,dt2,wdt2[3],wdt4[3];

  double **refcrd,ep=0.3,criteria=6.5;
  int NCmode=3,nibnum=3,numnb,num14;
  double *EAAm_equ,*ECGm_equ,EZm_equ;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char **inputfilename,**refcrdfilename,**parmfilename,**parmfilenameCG,*TACCMfilename;

  char **outputfilenameAAbase,**trjfilenameAAbase,**outputfilenameAA,**trjfilenameAA;
  char **outputfilenameCGbase,**trjfilenameCGbase,**outputfilenameCG,**trjfilenameCG;
  char *trjfileZbase,**trjfileThetaAAbase,**trjfileThetaCGbase,*logfilename;
  char trjfilenameZ[2000],**trjfilenameThetaAA,**trjfilenameThetaCG,logf[1000];
  FILE **inputfile,*refcrdfile,*parmfile,*TACCMfile,*logfile;

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
    {"numAA",1,NULL,'N'}, 
    {"numCG",1,NULL,'M'}, 
    {"tau",1,NULL,'a'}, {"dt",1,NULL,'x'},
    {"TAA",1,NULL,'T'}, 
    {"TCG",1,NULL,'1'},
    {"TZ",1,NULL,'B'},
    {"mZ",1,NULL,'m'}, {"massX",1,NULL,'X'},
    {"ep",1,NULL,'p'}, {"cutoff",1,NULL,'c'},
    {"AR",1,NULL,'A'},
    {"equ",1,NULL,'E'}, {"intequ",1,NULL,'I'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  TAA=(* double)gcemalloc(sizeof(double)*1);
  TCG=(* double)gcemalloc(sizeof(double)*1);

  i=0;j=0;

  while((c=getopt_long(argc,argv,"hv:n:s:e:a:i:x:N:M:T:1:B:m:X:p:c:A:E:I:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 's':
      numstep=atoi(optarg); break;
    case 'i':
      interval=atoi(optarg);  break;
    case 'e':
      numEX=atoi(optarg); break;
    case 'n':
      numRE=atoi(optarg); break;
    case 'N':
      numAA=atoi(optarg);
      TAA=(* double)gcerealloc(TAA,sizeof(double)*numAA);
      break;
    case 'M':
      numCG=atoi(optarg);  
      TCG=(* double)gcerealloc(TCG,sizeof(double)*numCG);
      break;
    case 'T':
      if (i<=numAA) {
	T0AA[i]=atof(optarg); 
	++i;
      }
      break;
    case '1':
      if (j<=numCG) {
	T0CG[j]=atof(optarg);  
	++j;
      }
      break;
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

  FGdata=(struct AADataforREMD_Amber *)gcerealloc(sizeof(struct AADataforMSREMD_Amber )*numAA);
  CGdata=(struct CGDataforREMD_PROTEINS2008 *)gcerealloc(sizeof(struct CGDataforMSREMD_PROTEINS2008 )*numCG);

  refcrdfilename=(char **)gcemalloc(sizeof(char *)*numCG);
  for (i=0;i<numCG;++i) refcrdfilename[i]=(char *)gcemalloc(sizeof(char)*2000);

  parmfilename=(char **)gcemalloc(sizeof(char *)numAA);
  for (i=0;i<numAA;++i) parmfilename[i]=(char *)gcemalloc(sizeof(char)*2000);

  outputfilenameAAbase=(char **)gcemalloc(sizeof(char *)numAA);
  for (i=0;i<numAA;++i) outputfilenameAAbase[i]=(char *)gcemalloc(sizeof(char)*2000);

  trjfilenameAAbase=(char **)gcemalloc(sizeof(char *)numAA);
  for (i=0;i<numAA;++i) trjfilenameAAbase[i]=(char *)gcemalloc(sizeof(char)*2000);

  parmfilenameCG=(char **)gcemalloc(sizeof(char *)numCG);
  for (i=0;i<numCG;++i) parmfilename[i]=(char *)gcemalloc(sizeof(char)*2000);

  outputfilenameCGbase=(char **)gcemalloc(sizeof(char *)numCG);
  for (i=0;i<numCG;++i) outputfilenameCGbase[i]=(char *)gcemalloc(sizeof(char)*2000);

  trjfilenameCGbase=(char **)gcemalloc(sizeof(char *)numCG);
  for (i=0;i<numCG;++i) trjfilenameCGbase[i]=(char *)gcemalloc(sizeof(char)*2000);

  ap=(struct AmberParmL *)gcemalloc(sizeof(struct AmberParmL)*numAA);
  e=(struct potential *)gcemalloc(sizeof(struct potential)*numAA);
  f=(struct force *)gcemalloc(sizeof(struct force)*numAA);

  apCG=(struct AmberParmL *)gcemalloc(sizeof(struct AmberParmL)*numCG);

  KTAA=(double *)gcemalloc(sizeof(double )*numAA);
  KTCG=(double *)gcemalloc(sizeof(double )*numCG);

  EAAm_equ=(double *)gcemalloc(sizeof(double)*numAA);
  ECGm_equ=(double *)gcemalloc(sizeof(double)*numCG);

  if (argc < 4+3*numAA+4*numCG) {
    USAGE(progname);
    exit(1);
  }
  inputfilename        = *argv;
  for (i=0;i<numCG;++i) refcrdfilename[i]    = *++argv;			   
  TACCMfilename        = *++argv;			 
  for (i=0;i<numAA;++i)  parmfilename[i]     = *++argv;
  for (i=0;i<numCG;++i)  parmfilenameCG[i]     = *++argv;
  for (i=0;i<numAA;++i)  outputfilenameAAbase[i] = *++argv;
  for (i=0;i<numAA;++i)  trjfilenameAAbase[i]    = *++argv;
  for (i=0;i<numCG;++i)  outputfilenameCGbase[i]= *++argv;
  for (i=0;i<numCG;++i)  trjfilenameCGbase   = *++argv;
  trjfileZbase         = *++argv;
  //  for (i=0;i<numAA;++i)   trjfileThetaAAbase[i]   = *++argv;
  //  for (i=0;i<numCG;++i)   trjfileThetaCGbase[i]   = *++argv;
  logfilename          = *++argv;

  index_replicas=(int *)gcemalloc(sizeof(int)*numRE);
  index_parameters=(int *)gcemalloc(sizeof(int)*numRE);
  REMD_ini_purmutation_funcs(numRE);

  for (i=0;i<numAA;++i) {
    parmfile=efopen(parmfilename[i],"r");
    readParmtopLb(parmfile,&ap[i]);
    fclose(parmfile); 
  }

  for (i=0;i<numCG;++i) {
    parmfile=efopen(parmfilenameCG[i],"r");
    readParmtopLb(parmfile,&ap_CG[i]);
    fclose(parmfile); 
  }

  /**********************************************************************************/
  /* Cdata.numatom=AP.NATOM;							    */
  /* j=0;  for (i=0;i<Cdata.numatom;++i)  if (strncmp(AP.IGRAPH[i],"H",1)==0)  ++j; */
  /* Cdata.numheavyatom=Cdata.numatom-j;  Cdata.numres=AP.NRES;			    */
  /* 										    */
  /* TACCMfile=efopen(TACCMfilename,"r");					    */
  /* fscanf(TACCMfile,"%d",&Zdata.numZ);					    */
  /* Zdata.pairs=(int **)gcemalloc(sizeof(int *)*Zdata.numZ);			    */
  /* for (i=0;i<Zdata.numZ;++i) Zdata.pairs[i]=(int *)gcemalloc(sizeof(int)*5);	    */
  /* for (i=0;i<Zdata.numZ;++i) {						    */
  /*   for (j=0;j<4;++j) fscanf(TACCMfile,"%d",&(Zdata.pairs[i][j]));		    */
  /*   fscanf(TACCMfile,"%d",&(Zdata.pairs[i][j]));				    */
  /* }										    */
  /* fclose(TACCMfile);								    */
  /**********************************************************************************/

  /*****************************************************************/
  /* Cdata.mass=(double *)gcemalloc(sizeof(double)*Cdata.numatom); */
  /* for (i=0;i<Cdata.numatom;++i) Cdata.mass[i]=AP.AMASS[i];	   */
  /*****************************************************************/
  refcrd=(double **)gcemalloc(sizeof(double *)*numCG);
  for (i=0;i<numCG;++i) refcrd[i]=(double *)gcemalloc(sizeof(double)*ap_CG[i].NATOM*3);

  for (i=0;i<numCG;++i) {
    refcrdfile=efopen(refcrdfilename[i],"r");
    getline(&line,&len,refcrdfile);
    fscanf(refcrdfile,"%d",&d);
    for (j=0;j<ap_CG[j].NATOM;++j) for (k=0;k<3;++k) fscanf(refcrdfile,"%lf",&refcrd[i][j*3+k]);
    fclose(refcrdfile);
  }

  for (i=0;i<numAA;++i) {
    FGdata[i].crd=(double *)gcemalloc(sizeof(double)*ap[i].NATOM*3);
    FGdata[i].vel=(double *)gcemalloc(sizeof(double)*ap[i].NATOM*3);
  }
  for (i=0;i<numCG;++i) {
    CGdata[i].crd=(double *)gcemalloc(sizeof(double)*ap_CG[i].NATOM*3);
    CGdata[i].vel=(double *)gcemalloc(sizeof(double)*ap_CG[i].NATOM*3);
  }

  Zdata.Z=(double *)gcemalloc(sizeof(double)*Zdata.numZ);
  Zdata.velZ=(double *)gcemalloc(sizeof(double)*Zdata.numZ);

  KZAA=(double **)gcemalloc(sizeof(double *)*numAA);
  for (i=0;i<numAA;++i) KZAA[i]=(double *)gcemalloc(sizeof(double)*numRE);
  KZCG=(double **)gcemalloc(sizeof(double *)*numCG);
  for (i=0;i<numCG;++i) KZCG[i]=(double *)gcemalloc(sizeof(double)*numRE);
  inputfile=efopen(inputfilename,"r");
  /***********************************************************************************************************/
  /* CGAAREMDreadInputs_Amber_PROTEINS2008_Amber_hybrid_1FG2CG(inputfile,Cdata.numatom,numRE,my_rank,	     */
  /* 							    FGdata.crd,FGdata.vel,			     */
  /* 							    CGdata[0].crd,CGdata[0].vel,		     */
  /* 							    CGdata[1].crd,CGdata[1].vel,		     */
  /* 							    KZAA,KZCG[0],KZCG[1]);			     */
  /***********************************************************************************************************/
  fclose(inputfile);
  for (i=0;i<numAA;++i) Zdata.KZAA[i]=KZAA[i][my_rank];
  for (i=0;i<numCG;++i) Zdata.KZCG[i]=KZCG[i][my_rank];

  if ( vMode==OFF ) {
    for (i=0;i<numAA;++i) {
      MD_Generate_inivelo(FGdata[i].vel,ap[i].AMASS,ap[i].NATOM,k_B*T0AA[i]*UNITT);
      FGdata[i].zeta=0.0;      FGdata[i].V_zeta=0.0;
    }
    for (i=0;i<numCG;++i) {
      MD_Generate_inivelo(CGdata[i].vel,ap_CG[i].AMASS,ap_CG[i].NATOM,k_B*T0CG[0]*UNITT);
      for (j=0;j<ap_CG[i].numatom;++j)
	if (strncmp(ap_CG[i].IGRAPH[j],"H",1)==0)
	  for (k=0;k<3;++k)
	    CGdata[i].vel[j*3+k]=0.0;
      CGdata[i].zeta=0.0;      CGdata[i].V_zeta=0.0;
    }

    TACCM_MD_Generate_inivelo(Zdata.velZ,Zdata.massZ,Zdata.numZ,k_B*T0Z*UNITT);
    Zdata.zetaZ=0.0;      Zdata.V_zetaZ=0.0;
  }  
  //  TACCM_CTheta(FGdata.crd,Cdata.numatom,Zdata.Z,Zdata.numZ,Zdata.pairs,pi);

  tau=tau/2.0/pi;  tau2=tau*tau;  
  KBTZ=k_B*T0Z;  Zdata.NfKTZ=(Zdata.numZ+1)*KBTZ*UNITT;  Zdata.QZ=tau2*KBTZ*UNITT*Zdata.numZ;
  KTAA=(double *)gcemalloc(sizeof(double)*numAA);
  KTCG=(double *)gcemalloc(sizeof(double)*numCG);

  for (i=0;i<numAA;++i) {
    KTAA[i]=k_B*T0AA[i];  
    FGdata[i].NfKT=(3.0*ap[i].NATOM+1)*KTAA[i]*UNITT;  
    FGdata[i].Q=tau2*KTAA[i]*UNITT*(3.0*ap[i].NATOM);
  }
  
  for (i=0;i<numCG;++i) {
    KTCG[i]=k_B*T0CG[i];
    CGdata[i].NfKT=(3.0*ap_CG[i].numheavyatom+1)*KTCG[i]*UNITT;  
    CGdata[i].Q=tau2*KTCG[i]*UNITT*(3.0*ap_CG[i].numheavyatom);
  }

  for (i=0;i<numAA;++i) {
    ffLc_set_calcffandforce(&(FGdata[i].e),&(FGdata[i].f),ap[i]);
  }

  for (i=0;i<numCG;++i) {
    ffLc_set_calcffandforce(&e[i],&f[i],ap_CG[i]);
    ffLc_set_non_bonding_index_1(&numnb,&num14,ap_CG[i]);
    e[i].parm.numnb=numnb;
    e[i].parm.num14=num14;
    e[i].parm.indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
    e[i].parm.index14=(int *)gcemalloc(sizeof(int)*num14*2);
    ffLc_set_non_bonding_index_2(e[i].parm.indexnb,e[i].parm.index14,ap_CG[i]);

    GOLMAA_PROTEINS2008_ff_set_calcff_b(&(CGdata[i].e),refcrd[i],Cdata.numatom,Cdata.numres,
					/*FGdata.*/e.parm.indexnb,/*FGdata.*/e.parm.numnb,ep,nibnum,criteria);
  }

  MD_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);

  for (i=0;i<numAA;++i) {
    sprintf(outputfilenameAA[i],"%s_%d",outputfilenameAAbase[i],my_rank+1);
    sprintf(trjfilenameAA[i],"%s_%d",trjfilenameAAbase[i],my_rank+1);
  }

  for (i=0;i<numCG;++i) {
    sprintf(outputfilenameCG[i],"%s_%d",outputfilenameCGbase[i],my_rank+1);
    sprintf(trjfilenameCG[i],"%s_%d",trjfilenameCGbase[i],my_rank+1);
  }

  sprintf(trjfilenameZ,"%s_%d",trjfileZbase,my_rank+1);
  //  sprintf(trjfilenameThetaAA,"%s_%d",trjfileThetaAAbase,my_rank+1);
  //  sprintf(trjfilenameThetaCG1,"%s_%d",trjfileThetaCG1base,my_rank+1);

  sprintf(logf,"%s_%d_ex.log",logfilename,my_rank+1);

  for (i=0;i<numAA;++i) {
    myncL_create_def_MCD(trjfilenameAA[i],ap[i].NATOM,&(FGdata[i].nc_id_MCD));
    FGdata[i].outputfile=efopen(outputfilenameAA[i],"w");
  }

  for (i=0;i<numCG;++i) {
    myncL_create_def_MCD(trjfilenameCG[i],ap_CG[i].numatom,&(CGdata[0].nc_id_MCD));
    CGdata[0].outputfile=efopen(outputfilenameCG[i],"w");
  }


  Zdata.trjfileZ=efopen(trjfilenameZ,"w");
  Zdata.trjfilThetaAA=efopen(trjfilenameThetaAA,"w");
  Zdata.trjfilThetaCG1=efopen(trjfilenameThetaCG1,"w");
  Zdata.trjfilThetaCG2=efopen(trjfilenameThetaCG2,"w");

  MPI_Barrier(MPI_COMM_WORLD);

  if ( num_procs != numRE ) {    printf("condition error\n");    exit(1);  }

  logfile=efopen(logf,"w");

  if (equflag==ON) {
    run_multi_MS_TACCM_Amber_PROTEINS2008_pca(numAA,numCG,
					      FGdata,CGdata,Zdata,
					      ap,ap_CG,
					      numstepequ,intervalequ,&l,
					      dt,dt2,wdt2,wdt4,nc,UNITT,k_B,pi,
					      EAAm_equ,ECGm_equ,&EZm_equ);

    for (i=0;i<numAA;++i) {
      FGdata[i].zeta=0.0;
      FGdata[i].V_zeta=0.0;
    }
    for (i=0;i<numCG;++i) {
      CGdata[i].zeta=0.0;
      CGdata[i].V_zeta=0.0;
    }
    Zdata.zetaZ=0.0; Zdata.V_zetaZ=0.0;
  }

  acc_ratio=MPI_MSREM_MD_MP1998_Amber_PROTEINS2008_pca(my_rank,num_procs,tag,&status,
						       numRE,numEX,KZAA,KZCG,
						       FGdata,CGdata,Zdata,Cdata,
						       T0AA,T0CG,T0Z,
						       numstep,interval,
						       dt,dt2,wdt2,wdt4,nc,
						       UNITT,k_B,tau,pi,
						       logfile);
  
  fclose(logfile);
  for (i=0;i<numAA;++i) {
    fclose(FGdata[i].outputfile);  
  }
  for (i=0;i<numCG;++i) {
    fclose(CGdata[i].outputfile);  
  }

  fclose(Zdata.trjfileZ);  

  if ( my_rank==0 ) {
    acc_ratio_file=efopen(acc_ratio_filename,"w");
    for (j=0;j<numRE;++j) {
      if (j!=numRE-1) fprintf(acc_ratio_file,"%d-%d %8.4lf\n",j,j+1,acc_ratio[j][j+1]);
      else fprintf(acc_ratio_file,"%d-0 %8.4lf\n",j,acc_ratio[j][0]);
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
