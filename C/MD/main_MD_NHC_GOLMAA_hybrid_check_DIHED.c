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
#include "FFL.h"
#include "EF.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

double MD_GOLMAAff_hybrid_calcff_check_DIHED(double *crd, int numatom,
					     FILE *parmtop,
					     int numsdihed,
					     int atom[4],
					     double dx,
					     double *f_d);

int main(int argc, char *argv[]) {
  int i,j,k,l,m,d;
  int numatom;
  double dx=0.0001;
  int numsdihed=10;
  int atom1,atom2,atom3,atom4;

  int DOF;

  double *f_d,*f_d2,p_d;
  int atom[4];

  double *crd,*refcrd,*mass;
  struct potential e;
  struct force f;
  struct potential_GOLMAA_hybrid e_GOLM;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*refcrdfilename,*parmfilename;

  FILE *inputfile,*reffile,*parmfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"dx",1,NULL,'x'},
    {"nums",1,NULL,'s'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hx:s:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'x':
      dx=atof(optarg);
      break;
    case 's':
      numsdihed=atoi(optarg);
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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  refcrdfilename = *++argv;
  parmfilename   = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      fscanf(inputfile,"%lf",&crd[i*3+j]);
    }
  }
  fclose(inputfile);

  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);
  reffile=efopen(refcrdfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(reffile,"%d",&d);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      fscanf(reffile,"%lf",&refcrd[i*3+j]);
    }
  }
  fclose(reffile);

  ffL_set_calcffandforce(&e,&f);

  f_d=(double *)gcemalloc(sizeof(double)*4*3);
  ffL_calcffandforce_speD(f_d,&p_d,crd,numsdihed,&atom1,&atom2,&atom3,&atom4);

  atom[0]=atom1;
  atom[1]=atom2;
  atom[2]=atom3;
  atom[3]=atom4;

  f_d2=(double *)gcemalloc(sizeof(double)*4*3);
  parmfile=efopen(parmfilename,"r");
  MD_GOLMAAff_hybrid_calcff_check_DIHED(crd,numatom,
					parmfile,
					numsdihed,
					atom,
					dx,f_d2);
  fclose(parmfile);

  printf("p_d=%e \n",p_d);

  for (i=0;i<4;++i) {
    printf("f_d (%d) =( ",i+1);
    for (j=0;j<3;++j) {
      printf("%e ,",f_d[i*3+j]);
    }
    printf(")\n");

    printf("f_d (%d) =( ",i+1);
    for (j=0;j<3;++j) {
      printf("%e ,",f_d2[i*3+j]);
    }
    printf(")\n");
    printf("\n");
  }

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename refcrdfilename clustfilename parmfilename outputfilename\n",progname);
}


double MD_GOLMAAff_hybrid_calcff_check_DIHED(double *crd, int numatom,
					     FILE *parmtop,
					     int numsdihed,
					     int atom[4],
					     double dx,
					     double *f_d) {
  int i,j,n;
  int atom1,atom2,atom3,atom4;
  double *crddx,*crddy,*crddz;
  double p_d,p_d_dx,p_d_dy,p_d_dz;
  double *f_d_dummy;

 readParmtopL(parmtop);

 crddx=(double *)gcemalloc(sizeof(double)*numatom*3);
 crddy=(double *)gcemalloc(sizeof(double)*numatom*3);
 crddz=(double *)gcemalloc(sizeof(double)*numatom*3);

 f_d_dummy=(double *)gcemalloc(sizeof(double)*12);
 ffL_calcffandforce_speD(f_d_dummy,&p_d,crd,numsdihed,&atom1,&atom2,&atom3,&atom4);

  for (n=0;n<4;++n) {
    for (i=0;i<numatom;++i) {
      for (j=0;j<3;++j) {
	crddx[i*3+j]=crd[i*3+j];
	crddy[i*3+j]=crd[i*3+j];
	crddz[i*3+j]=crd[i*3+j];
      }
    }

    crddx[abs(atom[n])]+=dx;
    crddy[abs(atom[n])+1]+=dx;
    crddz[abs(atom[n])+2]+=dx;
    
    ffL_calcffandforce_speD(f_d_dummy,&p_d_dx,crddx,numsdihed,&atom1,&atom2,&atom3,&atom4);
    ffL_calcffandforce_speD(f_d_dummy,&p_d_dy,crddy,numsdihed,&atom1,&atom2,&atom3,&atom4);
    ffL_calcffandforce_speD(f_d_dummy,&p_d_dz,crddz,numsdihed,&atom1,&atom2,&atom3,&atom4);

    f_d[n*3+0]=-(p_d_dx-p_d)/dx*UNIT;
    f_d[n*3+1]=-(p_d_dy-p_d)/dx*UNIT;
    f_d[n*3+2]=-(p_d_dz-p_d)/dx*UNIT;
  }
}
