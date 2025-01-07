#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "GOLM_Clementi_set.h"
#include "GOLM_Clementi.h"

#include "PTL.h"
#include "EF.h"
#include "RAND.h"
#include "BOXMULL.h"
#include "MD.h"

#include "netcdf_mineL.h"

#define NVT 1
#define NVE 0

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d;
  int numatom,numCAatom,numres,numstep=10000,interval=100;
  double dt=0.001;
  double *frc,PE;

  int vMode=OFF,MODE=NVT;

  double T0=300,T,K0,KE;
  double k_B=1.98723e-3;
  double UNITT=418.4070;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;

  double *crd,*refcrd,*refcrdAA,*mass,*vel,*vel0;

  double summass,COM[3];

  struct potential e_AA;
  struct potential_GOLM_Clementi e;

  double x[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*reffilename,*velfilename,*parmfilename;
  char *trjfilename,*outputfilename,*outputfilename2,*rstfilename="rstcrd",*rstvelfilename="rstvel";

  FILE *inputfile,*reffile,*velfile,*parmfile;
  FILE *outputfile,*outputfile2,*rstfile,*rstvelfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"nve",0,NULL,'*'},
    {"vMode",1,NULL,'v'},
    {"nums",1,NULL,'s'},
    {"temp",1,NULL,'t'},
    {"int",1,NULL,'i'},
    {"rst",1,NULL,'{'},
    {"rstvel",1,NULL,'}'},
    {"dt",1,NULL,'x'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"*hs:v:t:i:x:{:}:",long_opt,&opt_idx))!=-1) {
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
  reffilename       = *++argv;
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
  refcrd=(double *)gcemalloc(sizeof(double)*numCAatom*3);
  refcrdAA=(double *)gcemalloc(sizeof(double)*numatom*3);
  vel=(double *)gcemalloc(sizeof(double)*numCAatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  j=0;
  for (i=0;i<numatom;++i) {
    for (k=0;k<3;++k) fscanf(inputfile,"%lf",&x[k]);
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      for (k=0;k<3;++k) crd[j*3+k]=x[k];
      ++j;
    }
  }
  fclose(inputfile);

  reffile=efopen(reffilename,"r");
  getline(&line,&len,reffile);
  fscanf(reffile,"%d",&d);
  j=0;
  for (i=0;i<numatom;++i) {
    for (k=0;k<3;++k) fscanf(reffile,"%lf",&refcrdAA[i*3+k]);
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      for (k=0;k<3;++k) refcrd[j*3+k]=refcrdAA[i*3+k];
      ++j;
    }
  }
  fclose(reffile);

  if ( vMode==OFF ) {
    MD_Generate_inivelo(vel,mass,numCAatom,k_B*T0*UNITT);
  }
  else {
    velfile=efopen(velfilename,"r");
    for (i=0;i<numCAatom;++i) for (j=0;j<3;++j) fscanf(velfile,"%lf",&vel[i*3+j]);
    fclose(velfile);
  }
  K0=0.0;
  for (i=0;i<numCAatom;++i) 
    for (j=0;j<3;++j) 
      K0+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  for (i=0;i>numCAatom*3;++i) vel0[i]=vel[i];
  T=K0/((3*numCAatom-1.0)*k_B)*2.0/UNITT;

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
  GOLM_Clementi_ff_set_calcff(&e,refcrd,refcrdAA,numCAatom,numatom);
  GOLM_Clementi_ff_calcff(crd,numCAatom,&e);

  myncL_create_def_MCD(trjfilename,numCAatom,&nc_id_MCD);
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<numstep;++i) {
    if (MODE==NVT)
      KE=MD_Propagetor_Iso_JCP2003_GOLM_Clementi(crd,vel,mass,numCAatom,K0,dt,&e);
    else
      KE=MD_Propagetor_vV_NVE_GOLM_Clementi(crd,vel,mass,numCAatom,dt,&e);

    if (i%interval==0) {
      KE=KE/UNITT;
      T=KE/((3*numCAatom-1.0)*k_B)*2.0;

      fprintf(outputfile,"%d %e %e %e %e\n",i+1,e.p_t,KE,T,e.p_natatt_t);
      fprintf(outputfile2,"E_t    = %e \n",e.p_t+KE);
      fprintf(outputfile2,"KE     = %e \n",KE);
      fprintf(outputfile2,"p_tot  = %e \n",e.p_t);
      fprintf(outputfile2,"p_nat  = %e \n",e.p_natatt_t);
      fprintf(outputfile2,"p_rep  = %e \n",e.p_repul_t);
      fprintf(outputfile2,"p_dih  = %e \n",e.p_d_t);
      fprintf(outputfile2,"p_ang  = %e \n",e.p_a_t);
      fprintf(outputfile2,"p_bon  = %e \n",e.p_b_t);
      fprintf(outputfile2,"T      = %e \n",T);

      for (j=0;j<3;++j) COM[j]=0.0;
      for (j=0;j<numCAatom;++j) 
	for (k=0;k<3;++k) 
	  COM[k]+=mass[j]*crd[j*3+k]/summass;
      for (j=0;j<numCAatom;++j) 
	for (k=0;k<3;++k) 
	  crd[j*3+k]-=COM[k];

      for (j=0;j<numCAatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crd[j*3+k];
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
