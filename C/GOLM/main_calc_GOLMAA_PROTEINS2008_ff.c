#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"
#include "GOLMAA_PROTEINS2008_check.h"

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

  double f_natatt[3],f_repul[3],f_d[3],f_a[3],f_b[3];
  double dx=0.00001;
  int numspatom=10,nums=3,numa=4;

  double f_d1[4][3],f_d2[4][3];
  double f_natatt1[2][3],f_natatt2[2][3];

  double f_as[3][3],p_as0,f_as0[3][3];

  int **ncmap;
  double *crd,*refcrd;

  double ep=ep_natatt_PROTEINS2008;

  int nibnum=3,criteria=criteria_PROTEINS2008;

  struct potential e;
  struct force f;
  struct potential_GOLMAA_PROTEINS2008 e_GOLM;

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
    {"nums",1,NULL,'s'},
    {"numa",1,NULL,'a'},
    {"ep",1,NULL,'e'},
    {"cutoff",1,NULL,'c'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hx:s:a:e:c:",long_opt,&opt_idx))!=-1) {
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  inputfilename     = *argv;
  reffilename       = *++argv;
  parmfilename      = *++argv;
  outputfilename    = *++argv;
  //  mapfilename       = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  numres=AP.NRES;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  reffile=efopen(reffilename,"r");
  getline(&line,&len,reffile);
  fscanf(reffile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(reffile,"%lf",&refcrd[i*3+j]);
  fclose(reffile);

  ffL_set_calcffandforce(&e,&f);
  GOLMAA_PROTEINS2008_ff_set_calcff(&e_GOLM,refcrd,numatom,numres,e.parm.indexnb,e.parm.numnb,ep,nibnum,criteria);

  GOLMAA_PROTEINS2008_ff_calcff(crd,numatom,&e_GOLM);
  GOLMAA_PROTEINS2008_calcff_check(inputfilename,reffilename,parmfilename,
				   numspatom,dx,ep,nibnum,criteria,
				   f_natatt,f_repul,f_d,f_a,f_b,f_d1,f_d2,nums,
				   f_as,numa);

  GOLMAA_PROTEINS2008_ff_calcANGLE_forcheck(crd,numatom,(e_GOLM).Ka,(e_GOLM).ang_equ,
					    (e_GOLM).pairs_angl,(e_GOLM).num_angl,
					    &p_as0,f_as0,numa);

  printf("p_t=%e \n",e_GOLM.p_t);
  printf("p_nat=%e \n",e_GOLM.p_natatt_t);
  printf("p_rep=%e \n",e_GOLM.p_repul_t);
  printf("p_dih=%e \n",e_GOLM.p_d_t);
  printf("p_ang=%e \n",e_GOLM.p_a_t);
  printf("p_bon=%e \n",e_GOLM.p_b_t);

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
    printf("%e ,",e_GOLM.f_d[numspatom][i]);
  }
  printf(")\n");

  printf("f_d =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_d[i]);
  }
  printf(")\n");

  printf("f_a =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",e_GOLM.f_a[numspatom][i]);
  }
  printf(")\n");

  printf("f_a =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_a[i]);
  }
  printf(")\n");

  printf("f_b =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",e_GOLM.f_b[numspatom][i]);
  }
  printf(")\n");

  printf("f_b =( ");
  for (i=0;i<3;++i) {
    printf("%e ,",f_b[i]);
  }
  printf(")\n");

  for (i=0;i<3;++i) {
    printf("f_a %d =( ",numa);
    for (j=0;j<3;++j) {
      printf("%e ,",f_as0[i][j]);
    }
    printf(")\n");

    printf("f_a %d =( ",numa);
    for (j=0;j<3;++j) {
      printf("%e ,",f_as[i][j]);
    }
    printf(")\n");
  }


  /*****************************************************************************************************/
  /* ncmap=GOLMAA_PROTEINS2008_make_native_contact(refcrd,6.5,&((e_GOLM).num_natatt),numatom,numatom); */
  /* mapfile=efopen(mapfilename,"w");								       */
  /* for (i=0;i<numatom;++i) {									       */
  /*   for (j=0;j<numatom;++j) {								       */
  /*     fprintf(mapfile,"%2d ",ncmap[i][j]);							       */
  /*   }											       */
  /* }												       */
  /* fclose(mapfile);										       */
  /*****************************************************************************************************/


  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfilename parmfilename outputfilename\n",progname);
}


