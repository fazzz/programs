#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PT.h"
#include "FF.h"
#include "EF.h"
#include "IOV2.h"
#include "EF.h"
#include "quaternion.h"

#include "netcdf_mine.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,d;
  int numstep,numatom;
  double dt=0.001;
  double U=418.4070;

  double *crd,**v,*m;
  double *ene,**frc,**f_pre,ke;
  struct potential e;
  struct force f;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*outputfilename,*parmfilename;
  FILE *inputfile,*outputfile,*parmfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"h",long_opt,&opt_idx))!=-1) {
    switch(c) {
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
  numstep = atoi(*argv);
  inputfilename = *++argv;
  parmfilename = *++argv;
  outputfilename = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
 
  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  v=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i)v[i]=(double *)gcemalloc(sizeof(double)*3);
  m=(double *)gcemalloc(sizeof(double)*numatom);
  frc=(double **)gcemalloc(sizeof(double*)*numatom);
  for (i=0;i<numatom;++i) frc[i]=(double *)gcemalloc(sizeof(double)*3);
  f_pre=(double **)gcemalloc(sizeof(double*)*numatom);
  for (i=0;i<numatom;++i)f_pre[i]=(double *)gcemalloc(sizeof(double)*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      v[i][j]=0.0;
      f_pre[i][j]=0.0;
    }
  }
  for (i=0;i<numatom;++i) m[i]=AP.AMASS[i];

  ff_set_calcffandforce(&e,&f);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    ff_calcffandforce(crd,numatom,&e,&f);
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k) {
	frc[j][k]=-f.f_b[j*3+k]+f.f_a[j*3+k]+f.f_d[j*3+k]
	  +f.f_e[j*3+k]+f.f_e_14[j*3+k]+f.f_LJ[j*3+k]+f.f_LJ_14[j*3+k];
	v[j][k]+=dt/m[j]*(frc[j][k]+f_pre[j][k])/2.0;
	crd[j*3+k]+=dt*v[j][k]+dt*dt/2.0/m[j]*(frc[j][k]);
      }
    }
    ke=0.0;
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k) {
	f_pre[j][k]=frc[j][k];
	ke+=0.5*m[j]*v[j][k]*v[j][k];
      }
    }
    ke=ke/U;
    fprintf(outputfile,"%8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf\n",e.p_t,e.p_e_t,e.p_e_14_t,e.p_LJ_t,e.p_LJ_14_t,e.p_d_t,e.p_a_t,e.p_b_t);
  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parmfilename outputfilename \n",progname);
}
