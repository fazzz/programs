#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "GOLMAA_hybrid_set.h"
#include "GOLMAA_hybrid.h"

#include "PTL.h"
#include "EF.h"
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

  double ep=ep_natatt_hybrid;

  int vMode=OFF,MODE=NVT,NCmode=3,nibnum=1,criteria=criteria_hybrid;

  int nc=1;
  double T0=300,T,K0,KE;
  double k_B=1.98723e-3;
  double UNITT=418.4070;
  double NfKT,KT;
  double zeta=0.0,V_zeta=0.0,Q,tau=0.01,tau2;
  double PEv,KEv;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;

  double *crd,*refcrd,*refcrdAA,*mass,*vel;

  double summass,COM[3];

  struct potential e;
  struct force f;
  struct potential_GOLMAA_hybrid e_GOLM;
  double p_t=0.0;

  double x[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*refcrdfilename,*velfilename,*parmfilename;
  char *trjfilename,*outputfilename,*outputfilename2,*rstfilename="rstcrd",*rstvelfilename="rstvel";

  FILE *inputfile,*refcrdfile,*velfile,*parmfile;
  FILE *outputfile,*outputfile2,*rstfile,*rstvelfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"nve",0,NULL,'*'},
    {"vMode",1,NULL,'v'},
    {"ep",1,NULL,'e'},
    {"nums",1,NULL,'s'},
    {"temp",1,NULL,'t'},
    {"tau",1,NULL,'a'},
    {"int",1,NULL,'i'},
    {"rst",1,NULL,'{'},
    {"rstvel",1,NULL,'}'},
    {"dt",1,NULL,'x'},
    {"nc",1,NULL,'N'},
    {"nibnum",1,NULL,'b'},
    {"cutoff",1,NULL,'c'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"*hs:v:e:t:a:i:x:N:b:c:{:}:",long_opt,&opt_idx))!=-1) {
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
    case 'e':
      ep=atof(optarg);
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
    case 'N':
      NCmode=atoi(optarg);
      break;
    case 'b':
      nibnum=atoi(optarg);
      break;
    case 'c':
      criteria=atof(optarg);
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

  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }
  inputfilename     = *argv;
  refcrdfilename    = *++argv;
  parmfilename      = *++argv;
  outputfilename    = *++argv;
  outputfilename2   = *++argv;
  trjfilename       = *++argv;

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
  numres=AP.NATOM;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];
  
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);
  vel=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  refcrdfile=efopen(refcrdfilename,"r");
  getline(&line,&len,refcrdfile);
  fscanf(refcrdfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(refcrdfile,"%lf",&refcrd[i*3+j]);
  fclose(refcrdfile);

  if ( vMode==OFF ) {
    MD_Generate_inivelo(vel,mass,numatom,k_B*T0*UNITT);
    for (i=0;i<numatom;++i) {
      if (strncmp(AP.IGRAPH[i],"H",1)==0) {
	for (j=0;j<3;++j) {
	  vel[i*3+j]=0.0;
	}
      }
    }
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
  T=K0/((3*numheavyatom)*k_B)*2.0/UNITT;

  pi=acos(-1.0);
  tau=tau/2.0/pi;         
  tau2=tau*tau;          
  KT=k_B*T0;
  NfKT=(3.0*numheavyatom+1)*KT*UNITT;
  Q=tau2*KT*UNITT*(3.0*numheavyatom);

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
  GOLMAA_hybrid_ff_set_calcff_6_wtune(&e_GOLM,refcrd,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);

  ffL_calcffandforce_14vdWDAB_woH(crd,numatom,&e,&f);
  GOLMAA_hyb_ff_calcff(crd,numatom,&e_GOLM);

  MD_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);
  myncL_create_def_MCD(trjfilename,numatom,&nc_id_MCD);
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<numstep;++i) {
    if (MODE==NVT)
      KE=MD_Propagetor_NH_MP1998_GOLMAA_JCTC2011_wonat(crd,vel,mass,&zeta,&V_zeta,Q,NfKT,numatom,&KEv,&PEv,dt,dt2,nc,wdt4,wdt2,&e,&f,&e_GOLM);
    else
      KE=MD_Propagetor_vV_NVE_GOLMAA_JCTC2011_wonat(crd,vel,mass,numatom,dt,&e,&f,&e_GOLM);

    if (i%interval==0) {
      KE=KE/UNITT;
      T=KE/((3*numheavyatom)*k_B)*2.0;
      PEv=PEv/UNITT;
      KEv=KEv/UNITT;

      p_t=e_GOLM.p_repul_t+0.5*e.p_LJ_14_t+e.p_d_t+e.p_a_t+e.p_b_t;
      fprintf(outputfile,"%d %e %e %e %e %e %e %e\n",i+1,p_t,KE,KEv,PEv,p_t+KE+KEv+PEv,T,e_GOLM.p_natatt_t);
      fprintf(outputfile2,"E_t    = %e \n",p_t+KE+KEv+PEv);
      fprintf(outputfile2,"KE     = %e \n",KE);
      fprintf(outputfile2,"KEv    = %e \n",KEv);
      fprintf(outputfile2,"PEv    = %e \n",PEv);
      fprintf(outputfile2,"p_tot  = %e \n",p_t);
      fprintf(outputfile2,"p_rep  = %e \n",e_GOLM.p_repul_t);
      fprintf(outputfile2,"p_14LJ = %e \n",0.5*e.p_LJ_14_t);
      fprintf(outputfile2,"p_dih  = %e \n",e.p_d_t);
      fprintf(outputfile2,"p_ang  = %e \n",e.p_a_t);
      fprintf(outputfile2,"p_bon  = %e \n",e.p_b_t);
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

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename refcrdfilename parmfilename outputfilename outputfilename2 trjfilename\n",progname);
}
