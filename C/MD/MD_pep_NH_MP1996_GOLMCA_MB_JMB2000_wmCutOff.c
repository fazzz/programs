#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "GOLM_Clementi_set_wmCutOff.h"
#include "GOLM_Clementi.h"
#include "GOLM_Clementi_MB_wmCutOff.h"

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
  int i,j,k,l=0,dummy;
  int numatom,numCAatom,numres,numstep=10000,interval=100;
  double dt=0.001,dt2,wdt2[3],wdt4[3];
  double *frc,PE;
  double pi;

  double ep=ep_natatt_Clementi;
  double de=1.0,d=1.0,d2;

  int vMode=OFF,MODE=NVT;

  int nc=1;
  double T0=300,T,K0,KE;
  double k_B=1.98723e-3;
  double UNITT=418.4070;
  double NfKT,KT;
  double zeta=0.0,V_zeta=0.0,Q,tau=0.01,tau2;
  double PEv,KEv;

  double cutoff=6.5;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;

  double *crd,*refcrd1,*refcrd2,*refcrdAA1,*refcrdAA2,*mass,*vel;

  double summass,COM[3];

  struct potential e_AA;
  struct potential_GOLM_Clementi_MB e;

  double x[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*reffilename1,*reffilename2,*velfilename,*parmfilename;
  char *trjfilename,*outputfilename,*outputfilename2,*rstfilename="rstcrd",*rstvelfilename="rstvel";

  FILE *inputfile,*reffile1,*reffile2,*velfile,*parmfile;
  FILE *outputfile,*outputfile2,*rstfile,*rstvelfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"nve",0,NULL,'*'},
    {"vMode",1,NULL,'v'},
    {"nums",1,NULL,'s'},
    {"temp",1,NULL,'t'},
    {"tau",1,NULL,'a'},
    {"ep",1,NULL,'e'},
    {"int",1,NULL,'i'},
    {"rst",1,NULL,'{'},
    {"rstvel",1,NULL,'}'},
    {"dt",1,NULL,'x'},
    {"de",1,NULL,'d'},
    //    {"d",1,NULL,'2'},
    {"dd",1,NULL,'2'},
    {"cutoff",1,NULL,'c'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"*hs:v:t:a:e:i:x:d:c:2:{:}:",long_opt,&opt_idx))!=-1) {
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
    case 'e':
      ep=atof(optarg);
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
    case 'd':
      d=atof(optarg);
      break;
    case '2':
      de=atof(optarg);
      break;
    case 'c':
      cutoff=atof(optarg);
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

  if (argc < 7) {
    USAGE(progname);
    exit(1);
  }
  inputfilename     = *argv;
  reffilename1      = *++argv;
  reffilename2      = *++argv;
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
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      ++j;
    }
  }
  numCAatom=j;
  numres=AP.NATOM;
  mass=(double *)gcemalloc(sizeof(double)*numCAatom);
  j=0;
  for (i=0;i<numatom;++i) {
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      mass[j]=AP.AMASS[i];
      ++j;
    }
  }
  
  crd=(double *)gcemalloc(sizeof(double)*numCAatom*3);
  refcrd1=(double *)gcemalloc(sizeof(double)*numCAatom*3);
  refcrd2=(double *)gcemalloc(sizeof(double)*numCAatom*3);
  refcrdAA1=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrdAA2=(double *)gcemalloc(sizeof(double)*numatom*3);
  vel=(double *)gcemalloc(sizeof(double)*numCAatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&dummy);
  j=0;
  for (i=0;i<numatom;++i) {
    for (k=0;k<3;++k) fscanf(inputfile,"%lf",&x[k]);
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      for (k=0;k<3;++k) crd[j*3+k]=x[k];
      ++j;
    }
  }
  fclose(inputfile);

  reffile1=efopen(reffilename1,"r");
  getline(&line,&len,reffile1);
  fscanf(reffile1,"%d",&dummy);
  j=0;
  for (i=0;i<numatom;++i) {
    for (k=0;k<3;++k) fscanf(reffile1,"%lf",&refcrdAA1[i*3+k]);
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      for (k=0;k<3;++k) refcrd1[j*3+k]=refcrdAA1[i*3+k];
      ++j;
    }
  }
  fclose(reffile1);

  reffile2=efopen(reffilename2,"r");
  getline(&line,&len,reffile2);
  fscanf(reffile2,"%d",&dummy);
  j=0;
  for (i=0;i<numatom;++i) {
    for (k=0;k<3;++k) fscanf(reffile2,"%lf",&refcrdAA2[i*3+k]);
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      for (k=0;k<3;++k) refcrd2[j*3+k]=refcrdAA2[i*3+k];
      ++j;
    }
  }
  fclose(reffile2);

  if ( vMode==OFF ) {
    MD_Generate_inivelo(vel,mass,numCAatom,k_B*T0*UNITT);
    zeta=0.0;
    V_zeta=0.0;
  }
  else {
    velfile=efopen(velfilename,"r");
    for (i=0;i<numCAatom;++i) for (j=0;j<3;++j) fscanf(velfile,"%lf",&vel[i*3+j]);
    if (MODE==NVT) fscanf(velfile,"%lf %lf",&zeta,&V_zeta);
    fclose(velfile);
  }
  K0=0.0;
  for (i=0;i<numCAatom;++i) 
    for (j=0;j<3;++j) 
      K0+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  T=K0/((3*numCAatom)*k_B)*2.0/UNITT;

  pi=acos(-1.0);
  tau=tau/2.0/pi;         
  tau2=tau*tau;          
  KT=k_B*T0;
  NfKT=(3.0*numCAatom+1)*KT*UNITT;
  Q=tau2*KT*UNITT*(3.0*numCAatom);

  summass=0.0;
  for (i=0;i<numCAatom;++i) summass+=mass[i];
  for (i=0;i<3;++i) COM[i]=0.0;
  for (i=0;i<numCAatom;++i) 
    for (j=0;j<3;++j) 
      COM[j]+=mass[i]*crd[i*3+j]/summass;
  for (i=0;i<numCAatom;++i) 
    for (j=0;j<3;++j) 
      crd[i*3+j]-=COM[j];

  frc=(double *)gcemalloc(sizeof(double)*numCAatom*3);
  //  GOLM_Clementi_ff_set_calcff2(&e,refcrd,refcrdAA,numCAatom,numatom,ep);
  GOLM_Clementi_MB_ff_set_calcff_wmCutOff(&e,refcrd1,refcrd2,refcrdAA1,refcrdAA2,numCAatom,numatom,ep,cutoff);

  d2=d*d;
  d2=d2*k_B;
  de=de*k_B;

  //  GOLM_Clementi_ff_calcff(crd,numCAatom,&e);
  GOLM_Clementi_MB_ff_calcff(crd,numCAatom,de,d2,&e);

  MD_Propagetor_NH_Single_set_MP1996(nc,dt,&dt2,wdt2,wdt4);
  //  myncL_create_def_MCD(trjfilename,numatom,&nc_id_MCD);
  myncL_create_def_MCD(trjfilename,numCAatom,&nc_id_MCD);
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<numstep;++i) {
    if (MODE==NVT)
      KE=MD_Propagetor_NH_MP1998_GOLM_Clementi_MB(crd,vel,mass,&zeta,&V_zeta,Q,NfKT,numCAatom,&KEv,&PEv,dt,dt2,nc,wdt4,wdt2,&e,de,d2);
    /*****************************************************************************/
    /* else									 */
    /*   KE=MD_Propagetor_vV_NVE_GOLM_Clementi_MB(crd,vel,mass,numCAatom,dt,&e); */
    /*****************************************************************************/

    if (i%interval==0) {
      KE=KE/UNITT;
      T=KE/((3*numCAatom)*k_B)*2.0;
      PEv=PEv/UNITT;
      KEv=KEv/UNITT;

      fprintf(outputfile,"%d %e %e %e %e %e %e \n",i+1,e.p_MB,KE,KEv,PEv,e.p_MB+KE+KEv+PEv,T);
      fprintf(outputfile2,"E_t    = %e \n",e.p_MB+KE+KEv+PEv);
      fprintf(outputfile2,"KE     = %e \n",KE);
      fprintf(outputfile2,"KEv    = %e \n",KEv);
      fprintf(outputfile2,"PEv    = %e \n",PEv);
      fprintf(outputfile2,"T      = %e \n",T);

      for (j=0;j<3;++j) COM[j]=0.0;
      for (j=0;j<numCAatom;++j) 
	for (k=0;k<3;++k) 
	  COM[k]+=mass[j]*crd[j*3+k]/summass;
      for (j=0;j<numCAatom;++j) 
	for (k=0;k<3;++k) 
	  crd[j*3+k]-=COM[k];

      /*********************************************************/
      /* k=0;						       */
      /* for (j=0;j<numatom;++j) {			       */
      /* 	if (strncmp(AP.IGRAPH[j],"CA",2)==0) {	       */
      /* 	  for (l=0;l<3;++l) crd_nc[j][l]=crd[k*3+l];   */
      /* 	  ++k;					       */
      /* 	}					       */
      /* 	else { 					       */
      /* 	  for (l=0;l<3;++l) crd_nc[j][l]=0.0;	       */
      /* 	}					       */
      /* }						       */
      /*********************************************************/

      for (j=0;j<numatom;++j)
	for (k=0;k<3;++k) crd_nc[j][k]=crd[j*3+k];

      myncL_put_crd_ene_MCD(nc_id_MCD,l,crd_nc,e_AA,0.0);
      ++l;
    }
  }
  fclose(outputfile);
  fclose(outputfile2);
  nc_close((nc_id_MCD.ncid));

  rstfile=efopen(rstfilename,"w");
  fprintf(rstfile,"ACE\n");
  fprintf(rstfile,"%d\n",numatom);
  j=0;
  for (i=0;i<numatom;++i) {
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      for (k=0;k<3;++k) fprintf(rstfile,"%e ",crd[j*3+k]);
      ++j;
      fprintf(rstfile,"\n");
    }
    else { 
      for (k=0;k<3;++k) fprintf(rstfile,"%e ",0.0);
      fprintf(rstfile,"\n");
    }
  }
  fclose(rstfile);

  rstvelfile=efopen(rstvelfilename,"w");
  for (i=0;i<numCAatom;++i) {
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
