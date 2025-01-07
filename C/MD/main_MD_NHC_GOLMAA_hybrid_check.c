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

double MD_GOLMAAff_hybrid_calcff_check(double *crd, int numatom,
				       FILE *parmtop,
				       int numspatom,double dx,
				       double f_b[3],double f_a[3],double f_d[3],
				       double f_e_14[3],double f_LJ_14[3]);

int main(int argc, char *argv[]) {
  int i,j,k,l,m,d;
  int numatom;
  double dx=0.0001;
  int numspatom=21;
  double *frc,pot;

  int DOF;

  double f_b[3],f_a[3],f_d[3],f_e_14[3],f_LJ_14[3];

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_AMBER nc_id;

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
      numspatom=atoi(optarg);
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

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  ffL_set_calcffandforce(&e,&f);

  //  GOLMAA_hybrid_ff_set_calcff(&e_GOLM,refcrd,numatom,numres,e.parm.indexnb,e.parm.numnb);

  //  ffL_calcffandforce_14DAB_woH(crd,numatom,&e,&f);
  ffL_calcffandforce_woH(crd,numatom,&e,&f);

  parmfile=efopen(parmfilename,"r");
  MD_GOLMAAff_hybrid_calcff_check(crd,numatom,
				  parmfile,
				  numspatom,dx,
				  f_b,f_a,f_d,f_e_14,f_LJ_14);
  fclose(parmfile);

  printf("p_b=%e \n",e.p_b_t);
  printf("p_a=%e \n",e.p_a_t);
  printf("p_d=%e \n",e.p_d_t);
  printf("p_e_14=%e \n",e.p_e_14_t);
  printf("p_LJ_14=%e \n",e.p_LJ_14_t);

  printf("f_b =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f.f_b[numspatom*3+i]);
  }
  printf(")\n");

  printf("f_b =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_b[i]);
  }
  printf(")\n");
  printf("\n");

  printf("f_a =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f.f_a[numspatom*3+i]);
  }
  printf(")\n");

  printf("f_a =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_a[i]);
  }
  printf(")\n");
  printf("\n");

  printf("f_d =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f.f_d[numspatom*3+i]);
  }
  printf(")\n");

  printf("f_d =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_d[i]);
  }
  printf(")\n");
  printf("\n");

  printf("f_e14 =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f.f_e_14[numspatom*3+i]);
  }
  printf(")\n");
  printf("f_e14 =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_e_14[numspatom*3+i]);
  }
  printf(")\n");
  printf("\n");

  printf("f_LJ_14 =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f.f_LJ_14[numspatom*3+i]);
  }
  printf(")\n");
  printf("f_LJ_14 =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_LJ_14[i]);
  }
  printf(")\n");
  printf("\n");

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename refcrdfilename clustfilename parmfilename outputfilename\n",progname);
}


double MD_GOLMAAff_hybrid_calcff_check(double *crd, int numatom,
				       FILE *parmtop,
				       int numspatom,double dx,
				       double f_b[3],double f_a[3],double f_d[3],
				       double f_e_14[3],double f_LJ_14[3]) {
  int i,j;
  double *crddx,*crddy,*crddz;
  struct potential e,e_dx,e_dy,e_dz;
  struct force f,f_dx,f_dy,f_dz;

  readParmtopL(parmtop);

  crddx=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddy=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddz=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      crddx[i*3+j]=crd[i*3+j];
      crddy[i*3+j]=crd[i*3+j];
      crddz[i*3+j]=crd[i*3+j];
    }
  }

  crddx[numspatom*3]+=dx;
  crddy[numspatom*3+1]+=dx;
  crddz[numspatom*3+2]+=dx;

  ffL_set_calcffandforce(&e,&f);
  ffL_set_calcffandforce(&e_dx,&f_dx);
  ffL_set_calcffandforce(&e_dy,&f_dy);
  ffL_set_calcffandforce(&e_dz,&f_dz);

  ffL_calcffandforce_woH(crd,numatom,&e,&f);
  ffL_calcffandforce_woH(crddx,numatom,&e_dx,&f_dx);
  ffL_calcffandforce_woH(crddy,numatom,&e_dy,&f_dy);
  ffL_calcffandforce_woH(crddz,numatom,&e_dz,&f_dz);

  f_b[0]=-(e_dx.p_b_t-e.p_b_t)/dx*UNIT;
  f_b[1]=-(e_dy.p_b_t-e.p_b_t)/dx*UNIT;
  f_b[2]=-(e_dz.p_b_t-e.p_b_t)/dx*UNIT;

  f_a[0]=-(e_dx.p_a_t-e.p_a_t)/dx*UNIT;
  f_a[1]=-(e_dy.p_a_t-e.p_a_t)/dx*UNIT;
  f_a[2]=-(e_dz.p_a_t-e.p_a_t)/dx*UNIT;

  f_d[0]=-(e_dx.p_d_t-e.p_d_t)/dx*UNIT;
  f_d[1]=-(e_dy.p_d_t-e.p_d_t)/dx*UNIT;
  f_d[2]=-(e_dz.p_d_t-e.p_d_t)/dx*UNIT;

  f_e_14[0]=-(e_dx.p_e_14_t-e.p_e_14_t)/dx*UNIT;
  f_e_14[1]=-(e_dy.p_e_14_t-e.p_e_14_t)/dx*UNIT;
  f_e_14[2]=-(e_dz.p_e_14_t-e.p_e_14_t)/dx*UNIT;

  f_LJ_14[0]=-(e_dx.p_LJ_14_t-e.p_LJ_14_t)/dx*UNIT;
  f_LJ_14[1]=-(e_dy.p_LJ_14_t-e.p_LJ_14_t)/dx*UNIT;
  f_LJ_14[2]=-(e_dz.p_LJ_14_t-e.p_LJ_14_t)/dx*UNIT;

}
