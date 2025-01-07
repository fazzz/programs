#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "GOLMAA_hybrid_set.h"
#include "GOLMAA_hybrid.h"
#include "GOLMAA_hybrid_check2.h"
#include "FFL.h"

#include "PTL.h"
#include "EF.h"
#include "NC.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numatom,numres;

  double ep=ep_natatt_hybrid;

  int vMode=OFF,NCmode=3,nibnum=1,criteria=criteria_hybrid,CheckMode=OFF;

  double f_natatt[3],f_repul[3],f_d[3],f_a[3],f_b[3];
  double dx=0.00001;
  int numspatom=10,nums=3,numa=4;

  double f_d1[4][3],f_d2[4][3];
  double f_natatt1[2][3],f_natatt2[2][3];

  int **ncmap;
  double dt=0.001;
  double *crd,*refcrd,*mass,*frc;
  struct potential e;
  struct force f;
  struct potential_GOLMAA_hybrid e_GOLM;

  double x[3];

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*reffilename,*outputfilename,*parmfilename,*mapfilename;
  FILE *inputfile,*reffile,*outputfile,*parmfile,*mapfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"dx",1,NULL,'x'},
    {"ep",1,NULL,'e'},
    {"nums",1,NULL,'s'},
    {"numa",1,NULL,'a'},
    {"nc",1,NULL,'N'},
    {"nibnum",1,NULL,'b'},
    {"cutoff",1,NULL,'c'},
    {"check",1,NULL,'C'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hcx:e:s:a:N:b:c:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'x':
      dx=atof(optarg);
      break;
    case 's':
      numspatom=atoi(optarg);
      break;
    case 'a':
      numa=atoi(optarg);
      break;
    case 'e':
      ep=atof(optarg);
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
    case 'C':
      CheckMode=OFF;
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
  reffilename       = *++argv;
  parmfilename      = *++argv;
  outputfilename    = *++argv;
  mapfilename       = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  numres=AP.NRES;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  j=0;
  for (i=0;i<numatom;++i) mass[j]=AP.AMASS[i];

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  }
  fclose(inputfile);

  reffile=efopen(reffilename,"r");
  getline(&line,&len,reffile);
  fscanf(reffile,"%d",&d);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) fscanf(reffile,"%lf",&refcrd[i*3+j]);
  }
  fclose(reffile);

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  ffL_set_calcffandforce(&e,&f);
  //  GOLMAA_hybrid_ff_set_calcff(&e_GOLM,refcrd,numatom,numres,e.parm.indexnb,e.parm.numnb);

  //  ffL_calcffandforce_14vdWDAB_woH(crd,numatom,&e,&f);
  //  GOLMAA_hyb_ff_calcff(crd,numatom,&e_GOLM);
  GOLMAA_hybrid_ff_set_calcff_6_wtune(&e_GOLM,refcrd,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);

  ffL_calcffandforce_14vdWDAB_woH(crd,numatom,&e,&f);
  GOLMAA_hyb_ff_calcff(crd,numatom,&e_GOLM);

  printf("p_tot  = %e \n",e_GOLM.p_t+0.5*e.p_LJ_14_t+e.p_d_t+e.p_a_t+e.p_b_t);
  printf("p_nat  = %e \n",e_GOLM.p_natatt_t);
  printf("p_rep  = %e \n",e_GOLM.p_repul_t);
  printf("p_14LJ = %e \n",0.5*e.p_LJ_14_t);
  printf("p_dih  = %e \n",e.p_d_t);
  printf("p_ang  = %e \n",e.p_a_t);
  printf("p_bon  = %e \n",e.p_b_t);
  printf("\n");
  printf("p_nlo  = %e \n",e_GOLM.p_natatt_t);
  printf("p_loc  = %e \n",0.5*e.p_LJ_14_t+e.p_d_t+e.p_a_t+e.p_b_t);

  if (CheckMode==ON) {
    GOLM_AA_hybrid_calcff_check2(inputfilename,reffilename,parmfilename,numspatom,dx,
				 f_natatt,f_repul,f_d,f_a,f_b,nums);

    printf("p_t=%e \n",e_GOLM.p_t);
    printf("p_nat=%e \n",e_GOLM.p_natatt_t);
    printf("p_rep=%e \n",e_GOLM.p_repul_t);
    printf("p_dih=%e \n",e.p_d_t);
    printf("p_ang=%e \n",e.p_a_t);
    printf("p_bon=%e \n",e.p_b_t);

    printf("f_natt =( ");
    for (i=0;i<3;++i) {
      printf("%e ,",e_GOLM.f_natatt[numspatom][i]);
    }
    printf(")\n");

    printf("f_natt =( ");
    for (i=0;i<3;++i) {
      printf("%e ,",f_natatt[i]);
    }
    printf(")\n");

    printf("f_repl =( ");
    for (i=0;i<3;++i) {
      printf("%e ,",e_GOLM.f_repul[numspatom][i]);
    }
    printf(")\n");
    
    printf("f_repl =( ");
    for (i=0;i<3;++i) {
      printf("%e ,",f_repul[i]);
    }
    printf(")\n");

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
  }

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfilename parmfilename outputfilename\n",progname);
}


