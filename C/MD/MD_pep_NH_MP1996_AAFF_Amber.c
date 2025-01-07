#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

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
  double *frc,PE;
  double pi;

  int vMode=OFF,MODE=NVT,NCmode=3,nibnum=1;
  int EXCmode=OFF;
  int UMBdihedmode=OFF;

  int bodflag=ON,angflag=ON,dihflag=ON;
  int LJ14flag=ON,es14flag=ON;
  int LJflag=ON,esflag=ON;

  int nc=1;
  double T0=300,T,K0,KE;
  double k_B=1.98723e-3;
  double UNITT=418.4070;
  double NfKT,KT;
  double zeta=0.0,V_zeta=0.0,Q,tau=0.01,tau2;
  double PEv,KEv;

  int nump;

  double avePE=0.0,varPE=0.0,aveKE=0.0,varKE=0.0,aveT=0.0,varT=0.0;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;

  double *crd,*mass,*vel;

  double summass,COM[3];

  struct potential e;
  struct force f;
  double *pUMB,pUMB_t,**fUMB;
  double p_t=0.0,Etot;

  int *pairp;
  double *fcp,*dihe_equp;

  double x[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*velfilename,*parmfilename;
  char *trjfilename,*outputfilename,*outputfilename2,*rstfilename="rstcrd",*rstvelfilename="rstvel";
  char *umbfilename;

  char *logfilename="MD_pep_NH_MP1996_AAFF_Amber.log";

  FILE *inputfile,*velfile,*parmfile;
  FILE *outputfile,*outputfile2,*rstfile,*rstvelfile;
  FILE *umbfile;

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
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"*hbnf14LsEs:v:t:a:i:x:N:{:}:u:",long_opt,&opt_idx))!=-1) {
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
  inputfilename     = *argv;
  parmfilename      = *++argv;
  outputfilename    = *++argv;
  outputfilename2   = *++argv;
  trjfilename       = *++argv;

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
  
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  vel=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  if (UMBdihedmode==ON) {
    umbfile=efopen(umbfilename,"r");
    fscanf(umbfile,"%d",&nump);
    pairp=(int *)gcemalloc(sizeof(int)*nump*4);
    fcp=(double *)gcemalloc(sizeof(double)*nump);
    dihe_equp=(double *)gcemalloc(sizeof(double)*nump);
    pUMB=(int *)gcemalloc(sizeof(int)*nump);
    fUMB=(double **)gcemalloc(sizeof(double *)*numatom);
    for (i=0;i<numatom;++i)
      fUMB[i]=(double *)gcemalloc(sizeof(double)*3);
    for (i=0;i<nump;++i)
      fscanf(umbfile,"%d %d %d %d %lf %lf",&pairp[i*4+0],&pairp[i*4+1],&pairp[i*4+2],&pairp[i*4+3],&fcp[i],&dihe_equp[i]);
    fclose(umbfile);
    for (i=0;i<nump;++i) {
      for (j=0;j<4;++j) {
	pairp[i*4+j]-=1;
      }
      dihe_equp[i]=dihe_equp[i]*pi/180.0;
    }
  }

  if ( vMode==OFF ) {
    MD_Generate_inivelo(vel,mass,numatom,k_B*T0*UNITT);
    if (EXCmode==ON)
      for (i=0;i<numatom;++i)
	if (strncmp(AP.IGRAPH[i],"H",1)==0)
	  for (j=0;j<3;++j)
	    vel[i*3+j]=0.0;
    zeta=0.0;
    V_zeta=0.0;
  }
  else {
    velfile=efopen(velfilename,"r");
    for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(velfile,"%lf",&vel[i*3+j]);
    fclose(velfile);
  }
  K0=0.0;
  for (i=0;i<numatom;++i) 
    for (j=0;j<3;++j) 
      K0+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  if (EXCmode==OFF)
    T=K0/((3*numatom)*k_B)*2.0/UNITT;
  else
    T=K0/((3*numheavyatom)*k_B)*2.0/UNITT;

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

  summass=0.0;
  for (i=0;i<numatom;++i) summass+=mass[i];
  for (i=0;i<3;++i) COM[i]=0.0;
  for (i=0;i<numatom;++i) 
    for (j=0;j<3;++j) 
      COM[j]+=mass[i]*crd[i*3+j]/summass;
  for (i=0;i<numatom;++i) 
    for (j=0;j<3;++j) 
      crd[i*3+j]-=COM[j];

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  ffL_set_calcffandforce(&e,&f);

  /*if (EXCmode==OFF)*/ ffL_calcffandforce(crd,numatom,&e,&f);
  //  else ffL_calcffandforce_woH(crd,numatom,&e,&f);
  if (UMBdihedmode==ON) {
    pUMB_t=UMB_calc_dihetype_ff(crd,numatom,pairp,nump,fcp,dihe_equp,pUMB,fUMB);
  }

  MD_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);
  myncL_create_def_MCD(trjfilename,numatom,&nc_id_MCD);
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<numstep;++i) {
    if (MODE==NVT) {
      if (UMBdihedmode==OFF) {
	if (EXCmode==OFF)
	  KE=MD_Propagetor_NH_MP1998_AAFF_AMBER(crd,vel,mass,&zeta,&V_zeta,Q,NfKT,numatom,&KEv,&PEv,dt,dt2,nc,wdt4,wdt2,&e,&f);
	else
	  ;
	//	  KE=MD_Propagetor_NH_MP1998_AAFF_Amber_woH_wflag(crd,vel,mass,&zeta,&V_zeta,Q,NfKT,numatom,&KEv,&PEv,dt,dt2,nc,wdt4,wdt2,&e,&f,bodflag,angflag,dihflag,LJ14flag,es14flag,LJflag,esflag);
      }
      else {
	if (EXCmode==OFF)
	  KE=MD_Propagetor_NH_MP1998_AAFF_Amber_UMB(crd,vel,mass,&zeta,&V_zeta,Q,NfKT,numatom,&KEv,&PEv,dt,dt2,nc,wdt4,wdt2,&e,&f,pairp,nump,fcp,dihe_equp,pUMB,fUMB);
	else 
	  ;
	//	  KE=MD_Propagetor_NH_MP1998_AAFF_Amber_woH_UMB_wflag(crd,vel,mass,&zeta,&V_zeta,Q,NfKT,numatom,&KEv,&PEv,dt,dt2,nc,wdt4,wdt2,&e,&f,bodflag,angflag,dihflag,LJ14flag,es14flag,LJflag,esflag,pairp,nump,fcp,dihe_equp,fUMB);
      }
    }
    else {
      if (UMBdihedmode==OFF) {
	if (EXCmode==OFF)
	  KE=MD_Propagetor_vV_NVE_wflag(crd,vel,mass,numatom,dt,&e,&f,bodflag,angflag,dihflag,LJ14flag,es14flag,LJflag,esflag);
	else
	  ;
	//	KE=MD_Propagetor_vV_NVE_woH_wflag(crd,vel,mass,numatom,dt,&e,&f,bodflag,angflag,dihflag,LJ14flag,es14flag,LJflag,esflag);
      }
      else 
	if (EXCmode==OFF)
	  KE=MD_Propagetor_vV_NVE_UMB(crd,vel,mass,numatom,dt,&e,&f,pairp,nump,fcp,dihe_equp,pUMB,fUMB);
	else
	  ;
	//	KE=MD_Propagetor_vV_NVE_woH_wflag(crd,vel,mass,numatom,dt,&e,&f,bodflag,angflag,dihflag,LJ14flag,es14flag,LJflag,esflag);
    }
      
    if (i%interval==0) {
      KE=KE/UNITT;
      if (EXCmode==OFF)
	T=KE/((3*numatom)*k_B)*2.0;
      else
	T=KE/((3*numheavyatom)*k_B)*2.0;
      PEv=PEv/UNITT;
      KEv=KEv/UNITT;

      p_t=0.0;
      /* if (esflag==ON)    */p_t+=0.5*e.p_e_t;
      /* if (LJflag==ON)    */p_t+=0.5*e.p_LJ_t;
      /* if (es14flag==ON)	 */p_t+=0.5*e.p_e_14_t;
      /* if (LJ14flag==ON)	 */p_t+=0.5*e.p_LJ_14_t;
      /* if (dihflag==ON)	 */p_t+=e.p_d_t;
      /* if (angflag==ON)	 */p_t+=e.p_a_t;
      /* if (bodflag==ON)	 */p_t+=e.p_b_t;
      if (UMBdihedmode==ON) {
	pUMB_t=0.0; for (j=0;j<nump;++j) pUMB_t+=pUMB[j];
	p_t+=pUMB_t;
      }

      Etot=p_t+KE+KEv+PEv;
      fprintf(outputfile,"%d %e %e %e %e %e %e %e\n",i+1,p_t,KE,KEv,PEv,Etot,T);
      fprintf(outputfile2,"E_t    = %e \n",Etot);
      fprintf(outputfile2,"KE     = %e \n",KE);
      fprintf(outputfile2,"KEv    = %e \n",KEv);
      fprintf(outputfile2,"PEv    = %e \n",PEv);
      fprintf(outputfile2,"p_tot  = %e \n",/*e.*/p_t);
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
      ++l;
    }
  }
  fclose(outputfile);
  fclose(outputfile2);
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
