#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PTL.h"
#include "FFL.h"
#include "EF.h"
#include "MD.h"
#include "stringOTF.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

#define NVT 1
#define NVE 0

void usage(char *progname);
double Propagetor_Iso(double *crd,double *vel,int numatom,double dt,double *KE,double *PE);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,r,d;
  int numatom,numstep=10000,interval=100;
  int numcv=2;
  double dt=0.001,ang;
  double *frc,PE;
  double KE;
  double pi;

  int MODE=NVT;
  double Tobj=10,KEobj;
  double k_B=1.98723e-3,IsoCoff;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;

  double ***crd,***vel,*mass;
  struct potential e;
  struct force f;

  int numreplica=10;
  double **z,**z_star,*theta,**M;
  double kappa=1000,gamma=500,***force_eff;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *initialstringfname,*initialreplicfname,*parmfilename,*ouputstringfname;
  FILE *initialstring,*initialreplic,*parmfile,*ouputstring;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"nve",0,NULL,'*'},
    {"h",0,NULL,'h'},
    {"temp",1,NULL,'t'},
    {"t",1,NULL,'t'},
    {"int",1,NULL,'i'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"f*ht:i:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case '*':
      MODE=NVE;
      break;
    case 't':
      Tobj=atof(optarg);
      break;
    case 'i':
      interval=atoi(optarg);
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  initialstringfname = *argv;
  initialreplicfname = *++argv;
  parmfilename       = *++argv;
  ouputstringfname   = *++argv;

  pi=acos(-1.0);

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];

  crd=(double ***)gcemalloc(sizeof(double **)*2);
  vel=(double ***)gcemalloc(sizeof(double **)*2);
  for (i=0;i<2;++i) {
    crd[i]=(double **)gcemalloc(sizeof(double *)*numreplica);
    vel[i]=(double **)gcemalloc(sizeof(double *)*numreplica);
  }
  for (i=0;i<2;++i) {
    for (j=0;j<numreplica;++j) {
      crd[i][j]=(double *)gcemalloc(sizeof(double)*numatom*3);
      vel[i][j]=(double *)gcemalloc(sizeof(double)*numatom*3);
    }
  }
  z=(double **)gcemalloc(sizeof(double *)*numreplica);
  for (i=0;i<numreplica;++i) z[i]=(double *)gcemalloc(sizeof(double)*numcv);
  z_star=(double **)gcemalloc(sizeof(double *)*numreplica);
  for (i=0;i<numreplica;++i) z_star[i]=(double *)gcemalloc(sizeof(double)*numcv);
  theta=(double *)gcemalloc(sizeof(double)*numcv);
  M=(double **)gcemalloc(sizeof(double *)*numcv);
  for (i=0;i<numcv;++i)  M[i]=(double *)gcemalloc(sizeof(double)*numcv);

  initialstring=efopen(initialstringfname,"r");
  for (r=0;r<numreplica;++r) {
    for (i=0;i<numcv;++i) {
      fscanf(initialstring,"%lf",&ang);
      z[r][i]=ang*pi/180;
    }
  }
  fclose(initialstring);

  initialreplic=efopen(initialreplicfname,"r");
  for (i=0;i<numreplica;++i) {
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k) {
	fscanf(initialreplic,"%lf",&crd[0][i][j*3+k]);
	crd[1][i][j*3+k]=crd[0][i][j*3+k];
      }
    }
  }
  fclose(initialreplic);

  KEobj=0.5*3*numatom*k_B*Tobj/UNIT;
  IsoCoff=sqrt((3*numatom-1)*k_B*Tobj/UNIT);
  for (i=0;i<2;++i) for (j=0;j<numreplica;++j) MD_Generate_inivelo(vel[i][j],mass,numatom,k_B*Tobj/UNIT);

  ffL_set_calcffandforce(&e,&f);
  //  myncL_create_def_MCD(trjfilename,numatom,&nc_id_MCD);
  ouputstring=efopen(ouputstringfname,"w");
  for (i=0;i<numstep;++i) {
    for (r=0;r<numreplica;++r) {
      String_cTheta_FASYS(theta,crd[1][r]);
      String_cM_FASYS(crd[0][r],mass,M,numatom,numcv);
      String_Propagetor(z_star[r],z[r],numcv,M,theta,kappa,gamma,dt);

      String_Propagetor_Iso_FASYS(crd[0][r],vel[0][r],mass,numatom,IsoCoff,dt,&KE,&PE,e,f,z[r],kappa,numcv);
      String_Propagetor_Iso_FASYS(crd[1][r],vel[0][r],mass,numatom,IsoCoff,dt,&KE,&PE,e,f,z[r],kappa,numcv);
    }
    String_Repara(z,z_star,numreplica,numcv);

    if (i%interval==0) {
      for (j=0;j<numreplica;++j) {
	for (k=0;k<numcv;++k) fprintf(ouputstring,"%lf ",z[j][k]);
	fprintf(ouputstring,"\n");
      }
    }
  }
  fclose(ouputstringfname);
  //  nc_close((nc_id_MCD.ncid));

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parmfilename outputfilename\n",progname);
}

