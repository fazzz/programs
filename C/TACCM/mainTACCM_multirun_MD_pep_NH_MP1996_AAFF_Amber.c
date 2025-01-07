#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "TACCM.h"

#include "PTL.h"
#include "EF.h"
#include "UMBP.h"
#include "RAND.h"
#include "BOXMULL.h"
#include "MD.h"
#include "MD_NHC_MP1996.h"

#include "netcdf_mineL.h"

#define NVT 1
#define NVE 0

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d;
  int numatom,numheavyatom,numres,numstep=10000,interval=100;
  double dt=0.001,dt2,wdt2[3],wdt4[3];
  double **frc,*PE;  // Multi Run
  double pi;

  int vMode=OFF,MODE=NVT,NCmode=3,nibnum=1;
  int EXCmode=OFF;
  int UMBdihedmode=OFF;

  int bodflag=ON,angflag=ON,dihflag=ON;
  int LJ14flag=ON,es14flag=ON;
  int LJflag=ON,esflag=ON;
  ///////////////// TACCM //////////////////////
  int massflag=OFF;
  double massX=1.0;
  ///////////////// TACCM //////////////////////

  int nc=1;
  double T0=300,*T,*K0,*KE;             // Multi Run
  double k_B=1.98723e-3;
  double UNITT=418.4070;
  double NfKT,KT;
  double *zeta,*V_zeta,Q,tau=0.01,tau2; // Multi Run
  double *PEv,*KEv;                     // Multi Run

  int nump;

  double avePE=0.0,varPE=0.0,aveKE=0.0,varKE=0.0,aveT=0.0,varT=0.0;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD *nc_id_MCD; // Multi Run
  //  struct my_netcdf_out_id_AMBER nc_id_MCD;

  int numrun=1; // Multi Run

  double **crd,*mass,**vel; // Multi Run

  double summass,*COM[3]; // Multi Run

  struct potential *e; // Multi Run
  struct force *f; // Multi Run
  //  double *pUMB,pUMB_t,**fUMB;
  double *p_t=0.0,*Etot;

  int *pairp;
  double *fcp,*dihe_equp; // Multi Run


  ///////////////// TACCM //////////////////////
  double *theta;
  double *Z,*velZ,*accZ;
  int numZ;
  double TobjZ=500,KEobjZ,KBTZ,TZ;
  double massZ=100.0;
  double *frcZ;
  double KEZ,PEZ,KEvZ,PEvZ,EtZ;
  double zetaZ=0.0,V_zetaZ=0.0,QZ,NfKTZ;
  double KZ=10.0;
  int *indexTACCM,**pairsZ;
  char *TACCMfilename,*trjfilenameZ,*trjfilenameTheta;
  FILE *TACCMfile,*trjfileZ,*trjfileTheta;

  double tauZ,tau2Z;
  ///////////////// TACCM //////////////////////

  double x[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilelistname,*velfilelistname;       // Multi Run
  char *inputfilename,*velfilename,*parmfilename; // Multi Run
  char *trjfilenamebase,*outputfilenamebase;      // Multi Run
  char *rstfilenamebase,*rstvelfilenamebase;      // Multi Run
  char *trjfilename[1000],*outputfilename[1000];  // Multi Run
  char *rstfilename[1000],*rstvelfilename[1000];  // Multi Run
  //  char *umbfilename;

  char *logfilename="TACCM_nultirun_MD_pep_NH_MP1996_AAFF_Amber.log";

  FILE *inputfilelist,*velfilelist,*parmfile; // Multi Run
  FILE *inputfile,*velfile;                   // Multi Run
  FILE **outputfile,**rstfile,**rstvelfile;   // Multi Run
  //  FILE *umbfile;

  FILE *logfile;

  char *progname;

  int opt_idx=1;

  pi=acos(-1.0);

  struct option long_opt[] = {
    {"nve",0,NULL,'*'},
    {"vMode",1,NULL,'v'},
    {"nums",1,NULL,'s'},
    {"temp",1,NULL,'t'},
    {"tau",1,NULL,'a'},
    {"int",1,NULL,'i'},
    {"rst",1,NULL,'{'},
    {"rstvel",1,NULL,'}'},
    {"dt",1,NULL,'x'},
    {"BodOff",0,NULL,'b'},
    {"AngOff",0,NULL,'n'},
    {"DihOff",0,NULL,'f'},
    {"14LJOff",0,NULL,'1'},
    {"14esOff",0,NULL,'4'},
    {"LJOff",0,NULL,'L'},
    {"esOff",0,NULL,'E'},
    {"Hexc",0,NULL,'e'},
    {"UMB",1,NULL,'u'},
    // TACCM ///////////////
    {"mZ",1,NULL,'m'},
    {"KZ",1,NULL,'K'},
    {"tempB",1,NULL,'B'},
    {"massX",1,NULL,'X'},
    // TACCM ///////////////
    // Multi Run ///////////
    {"numrun",1,NULL,'M'}
    // Multi Run ///////////
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"*hbnf14LsEs:v:t:a:i:x:N:{:}:u:m:K:B:X:M:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case '*':
      MODE=NVE;
      break;
    case 's':
      numstep=atoi(optarg);
      break;
    case 'v':
      vMode=ON;
      velfilename=optarg;
      break;
    case 't':
      T0=atof(optarg);
      break;
    case 'a':
      tau=atof(optarg);
      break;
    case 'i':
      interval=atoi(optarg);
      break;
    case '{':
      rstfilename=optarg;
      break;
    case '}':
      rstvelfilename=optarg;
      break;
    case 'x':
      dt=atof(optarg);
      break;
    case 'b':
      bodflag=OFF;
      break;
    case 'n':
      angflag=OFF;
      break;
    case 'f':
      dihflag=OFF;
      break;
    case '1':
      LJ14flag=OFF;
      break;
    case '4':
      es14flag=OFF;
      break;
    case 'L':
      LJflag=OFF;
      break;
    case 'E':
      esflag=OFF;
      break;
    case 'e':
      EXCmode=ON;
      break;
    case 'u':
      umbfilename=optarg;
      UMBdihedmode=ON;
      break;
      ///////////////// TACCM //////////////////////
    case 'm':
      massZ=atof(optarg);
      break;
    case 'K':
      KZ=atof(optarg);
      break;
    case 'B':
      TobjZ=atof(optarg);
      break;
    case 'X':
      massflag=ON;
      massX=atof(optarg);
      break;
      ///////////////// TACCM //////////////////////
      // Multi Run ///////////
    case 'M':
      numrun=atoi(optarg);
      break;
    // Multi Run ///////////
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
  parmfilename      = *++argv;
  ///////////////// TACCM //////////////////////
  TACCMfilename  = *++argv;
  ///////////////// TACCM //////////////////////
  outputfilenamebase    = *++argv;              // Multi Run
  trjfilenamebase       = *++argv;              // Multi Run
  trjfilenameZbase   = *++argv;                 // Multi Run
  trjfilenameThetabase   = *++argv;             // Multi Run

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  j=0;
  for (i=0;i<numatom;++i) if (strncmp(AP.IGRAPH[i],"H",1)==0)  ++j;
  numheavyatom=numatom-j;
  numres=AP.NRES;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];
  if ( massflag==ON ) {
    for (i=0;i<numatom;++i) mass[i]=massX;
  }
  
  crd=(double **)gcemalloc(sizeof(double *)*numrun); // Multi Run
  vel=(double **)gcemalloc(sizeof(double *)*numrun); // Multi Run

  for (i=0;i<numrun;++i) {                                // Multi Run
    crd[i]=(double *)gcemalloc(sizeof(double)*numatom*3); // Multi Run
    vel[i]=(double *)gcemalloc(sizeof(double)*numatom*3); // Multi Run
  }                                                       // Multi Run

  // Multi Run ///////////
  inputfilelist=efopen(inputfilelistname,"r");
  for (i=0;i<numrun;++i) {
    getline(&line,&len,inputfilelist);
    for (j=0;line[j]!=NULL;++j)
      inputfilename[j]=line[j];
    inputfilename[j]=NULL;
    inputfile=efopen(inputfilename,"r");
    getline(&line,&len,inputfile);
    fscanf(inputfile,"%d",&d);
    for (j=0;j<numatom;++j) for (k=0;k<3;++k) fscanf(inputfile,"%lf",&crd[i][j*3+k]);
    fclose(inputfile);
  }
  fclose(inputfilelist);
  // Multi Run ///////////

  // Multi Run ///////////
  zeta=(double *)gcemalloc(sizeof(double)*numrun);
  V_zeta=(double *)gcemalloc(sizeof(double)*numrun);
  T=(double *)gcemalloc(sizeof(double)*numrun);
  K0=(double *)gcemalloc(sizeof(double)*numrun);
  KE=(double *)gcemalloc(sizeof(double)*numrun);
  PEv=(double *)gcemalloc(sizeof(double)*numrun);
  KEv=(double *)gcemalloc(sizeof(double)*numrun);
  if ( vMode==OFF ) {
    for (i=0;i<numrun:++i) {
      MD_Generate_inivelo(vel[i],mass,numatom,k_B*T0*UNITT);
      if (EXCmode==ON)
	for (j=0;j<numatom;++j)
	if (strncmp(AP.IGRAPH[j],"H",1)==0)
	  for (k=0;k<3;++k)
	    vel[i][j*3+k]=0.0;
      zeta[i]=0.0;
      V_zeta[i]=0.0;
    }
  }
  else {
    velfilelist=efopen(velfilelistname,"r");
    for (i=0;i<numrun:++i) {
      getline(&line,&len,velfilelist);
      for (j=0;line[j]!=NULL;++j)
	velfilename[j]=line[j];
      velfilename[j]=NULL;
      velfile=efopen(velfilename,"r");
      for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(velfile,"%lf",&vel[i][j*3+k]);
      fclose(velfile);
    }
  }
  fclose(velfilelistname);
  // Multi Run ///////////

  // Multi Run ///////////
  for (i=0;i<numrun;++i) {
    K0[i]=0.0;
    for (j=0;j<numatom;++j) 
      for (k=0;k<3;++k) 
	K0[i]+=0.5*mass[j]*vel[i][j*3+k]*vel[i][j*3+k];
    if (EXCmode==OFF)
      T[i]=K0[i]/((3*numatom)*k_B)*2.0/UNITT;
    else
      T[i]=K0[i]/((3*numheavyatom)*k_B)*2.0/UNITT;
  }
  // Multi Run ///////////

  tau=tau/2.0/pi;         
  tau2=tau*tau;          
  KT=k_B*T0;
  if (EXCmode==OFF) {
    NfKT=(3.0*numatom+1)*KT*UNITT;
    Q=tau2*KT*UNITT*(3.0*numatom);
  }
  else {
    NfKT=(3.0*numheavyatom+1)*KT*UNITT;
    Q=tau2*KT*UNITT*(3.0*numheavyatom);
  }

  // Multi Run ///////////
  summass=0.0;
  for (i=0;i<numatom;++i) summass+=mass[i];
  for (i=0;i<numrun;++i)  for (j=0;j<3;++j) COM[i][j]=0.0;
  for (i=0;i<numrun;++i)
    for (j=0;j<numatom;++j) 
      for (k=0;k<3;++k) 
	COM[i][j]+=mass[i]*crd[i][j*3+k]/summass;
  for (i=0;i<numrun;++i)
    for (j=0;j<numatom;++j) 
      for (k=0;k<3;++k) 
	crd[i][j*3+k]-=COM[i][k];
  // Multi Run ///////////

  frc=(double **)gcemalloc(sizeof(double *)*numrun);      // Multi Run
  for (i=0;i<numrun;++i)                                  // Multi Run
    frc[i]=(double *)gcemalloc(sizeof(double)*numatom*3); // Multi Run

  for (i=0;i<numrun;++i) ffL_set_calcffandforce(&(e[i]),&(f[i])); // Multi Run

  for (i=0;i<numrun;++i) ffL_calcffandforce(crd[i],numatom,&(e[i]),&(f[i]));   // Multi Run

  ///////////////// TACCM //////////////////////
  TACCMfile=efopen(TACCMfilename,"r");
  fscanf(TACCMfile,"%d",&numZ);
  pairsZ=(int **)gcemalloc(sizeof(int *)*numZ);
  for (i=0;i<numZ;++i) pairsZ[i]=(int *)gcemalloc(sizeof(int)*5);
  for (i=0;i<numZ;++i) {
    for (j=0;j<4;++j) 
      fscanf(TACCMfile,"%d",&pairsZ[i][j]);
    fscanf(TACCMfile,"%d",&pairsZ[i][j]);
  }
  fclose(TACCMfile);
  theta=(double *)gcemalloc(sizeof(double)*numZ);
  Z=(double *)gcemalloc(sizeof(double)*numZ);
  velZ=(double *)gcemalloc(sizeof(double)*numZ);
  frcZ=(double *)gcemalloc(sizeof(double)*numZ);
  TACCM_CTheta(crd,numatom,theta,numZ,pairsZ,pi);
  for (i=0;i<numZ;++i) Z[i]=theta[i];
  KBTZ=k_B*TobjZ;
  NfKTZ=(numZ+1)*KBTZ*UNITT;
  QZ=tau2*KBTZ*UNITT*numZ;
  KEZ=TACCM_MD_Generate_inivelo(velZ,massZ,numZ,k_B*TobjZ*UNITT);
  ///////////////// TACCM //////////////////////

  MD_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);

  nc_id_MCD=(struct my_netcdf_out_id_MCD *)gcemalloc(sizeof(struct my_netcdf_out_id_MCD)*numrun);
  trjfilename=(char )
  for (i=0;i<numrun;++i) {                                 // Multi Run
    myncL_create_def_MCD(trjfilename[i],numatom,&(nc_id_MCD[i]));
    //  myncL_create_def_AMBER(trjfilename,numatom,&nc_id_MCD);
    outputfile=efopen(outputfilename,"w");
    outputfile2=efopen(outputfilename2,"w");
    ///////////////// TACCM //////////////////////
    trjfileZ=efopen(trjfilenameZ,"w");
    trjfileTheta=efopen(trjfilenameTheta,"w");
    ///////////////// TACCM //////////////////////
  }

  for (i=0;i<numstep;++i) {
    TACCM_CTheta(crd,numatom,theta,numZ,pairsZ,pi);

    TACCM_MD_Propagetor_NH_MP1998_AAFF_Amber(crd,vel,mass,&zeta,&V_zeta,Q,NfKT,numatom,&KE,&KEv,&PEv,dt,dt2,nc,wdt4,wdt2,&e,&f,Z,numZ,theta,KZ,pairsZ,&PEZ,pi);

    TACCM_MD_Propagetor_NH_MP1998_Z(Z,velZ,massZ,theta,&zetaZ,&V_zetaZ,QZ,NfKTZ,numZ,&KEZ,&KEvZ,&PEvZ,dt,dt2,nc,wdt4,wdt2,KZ,&PEZ,frcZ,pi);
      
    if (i%interval==0) {
      KE=KE/UNITT;
      T=KE/((3*numatom)*k_B)*2.0;

      PEv=PEv/UNITT;
      KEv=KEv/UNITT;

      p_t=0.0;
      p_t+=0.5*e.p_e_t;
      p_t+=0.5*e.p_LJ_t;
      p_t+=0.5*e.p_e_14_t;
      p_t+=0.5*e.p_LJ_14_t;
      p_t+=e.p_d_t;
      p_t+=e.p_a_t;
      p_t+=e.p_b_t;

      Etot=p_t+KE+KEv+PEv+PEZ;
      fprintf(outputfile,"%d %e %e %e %e %e %e ",i+1,p_t,KE,KEv,PEv,Etot,T);
      fprintf(outputfile2,"E_t    = %e \n",Etot);
      fprintf(outputfile2,"KE     = %e \n",KE);
      fprintf(outputfile2,"KEv    = %e \n",KEv);
      fprintf(outputfile2,"PEv    = %e \n",PEv);
      fprintf(outputfile2,"p_tot  = %e \n",/*e.*/p_t);

      ///////////////// TACCM //////////////////////
      KEZ=KEZ/UNITT;
      TZ=KEZ/(numZ*k_B)*2.0;

      PEvZ=PEvZ/UNITT;
      KEvZ=KEvZ/UNITT;

      EtZ=PEZ+KEZ+PEvZ+KEvZ;
      fprintf(outputfile,"%d %e %e %e %e %e %e %e\n",i+1,PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);
      ///////////////// TACCM //////////////////////

      if (esflag==ON)	fprintf(outputfile2,"p_es   = %e \n",0.5*e.p_e_t);
      if (LJflag==ON)	fprintf(outputfile2,"p_LJ   = %e \n",0.5*e.p_LJ_t);
      if (es14flag==ON)	fprintf(outputfile2,"p_14es = %e \n",0.5*e.p_e_14_t);
      if (LJ14flag==ON)	fprintf(outputfile2,"p_14LJ = %e \n",0.5*e.p_LJ_14_t);
      if (dihflag==ON)	fprintf(outputfile2,"p_dih  = %e \n",e.p_d_t);
      if (angflag==ON)	fprintf(outputfile2,"p_ang  = %e \n",e.p_a_t);
      if (bodflag==ON)	fprintf(outputfile2,"p_bon  = %e \n",e.p_b_t);
      if (UMBdihedmode==ON) fprintf(outputfile2,"p_UMB  = %e \n",pUMB_t);
      fprintf(outputfile2,"T      = %e \n",T);

      for (j=0;j<3;++j) COM[j]=0.0;
      for (j=0;j<numatom;++j) 
	for (k=0;k<3;++k) 
	  COM[k]+=mass[j]*crd[j*3+k]/summass;
      for (j=0;j<numatom;++j) 
	for (k=0;k<3;++k) 
	  crd[j*3+k]-=COM[k];

      for (j=0;j<numatom;++j)
	for (k=0;k<3;++k) crd_nc[j][k]=crd[j*3+k];

      avePE=(i*avePE+e.p_t)/(i+1);
      varPE=(i*varPE+e.p_t*e.p_t)/(i+1);

      aveKE=(i*aveKE+KE)/(i+1);
      varKE=(i*varKE+KE*KE)/(i+1);

      aveT=(i*aveT+T)/(i+1);
      varT=(i*varT+T*T)/(i+1);

      myncL_put_crd_ene_MCD(nc_id_MCD,l,crd_nc,e,0.0);
      //      myncL_put_crd_AMBER(nc_id_MCD,l,crd_nc);
      ++l;

      ///////////////// TACCM //////////////////////
      for (j=0;j<numZ;++j) {
	fprintf(trjfileZ,"%e ",Z[j]);
      }
      fprintf(trjfileZ,"\n");
      for (j=0;j<numZ;++j) {
	fprintf(trjfileTheta,"%e ",theta[j]);
      }
      fprintf(trjfileTheta,"\n");      
      ///////////////// TACCM //////////////////////
    }
  }
  fclose(outputfile);
  fclose(outputfile2);
  ///////////////// TACCM //////////////////////
  fclose(trjfileZ);
  fclose(trjfileTheta);
  ///////////////// TACCM //////////////////////
  nc_close((nc_id_MCD.ncid));

  rstfile=efopen(rstfilename,"w");
  fprintf(rstfile,"ACE\n");
  fprintf(rstfile,"%d\n",numatom);
  for (i=0;i<numatom;++i) {
    for (k=0;k<3;++k) fprintf(rstfile,"%e ",crd[i*3+k]);
    fprintf(rstfile,"\n");
  }
  fclose(rstfile);

  rstvelfile=efopen(rstvelfilename,"w");
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) fprintf(rstvelfile,"%e ",vel[i*3+j]);
    if (MODE==NVT) fprintf(rstvelfile,"%lf %lf",zeta,V_zeta);
    fprintf(rstvelfile,"\n");
  }
  fclose(rstvelfile);

  logfile=efopen(logfilename,"w");
  varPE=sqrt(varPE);
  fprintf(logfile,"Potential energy = %10.5lf kcal/mol +- %10.5lf kcal/mol\n",avePE,varPE);
  varKE=sqrt(varKE);
  fprintf(logfile,"Kinetic   energy = %10.5lf kcal/mol +- %10.5lf kcal/mol\n",aveKE,varKE);
  varT=sqrt(varT);
  fprintf(logfile,"Temperature      = %10.5lf K        +- %10.5lf K\n",aveT,varT);
  fclose(logfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename refcrdfilename parmfilename outputfilename outputfilename2 trjfilename\n",progname);
}

