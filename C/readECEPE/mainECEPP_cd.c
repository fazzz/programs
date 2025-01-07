#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ECEPE.h"
#include "PT.h"
#include "FF.h"

#include "EF.h"

#define ON  1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int d;

  int numnb,num14;
  int *indexnb,*index14,*pairs;

  double *co;
  double *dihed,*delta_dihed;
  double pi;
  double dihed_check;

  double *co_dummy,*dihed_dummy,*ene;

  int *bp;
  int **bp_f;
  int **pair1_5,**pair1_4;
  int *num1_4,*num1_5;
  int *numb;
  int atom[4];

  struct ECEPE_parms ECEPE_p;
  struct pnb nb_p;
  struct ECEPE_pote p;
  struct ECEPE_force f;

  char *progname;
  char *preofilename,*bd8filename,*parmtopfilename,*indexfilename,*coofilename,*coofilename_for_sflag,*outfilename;
  FILE *parmtopfile,*indexfile,*coofile,*dp,*coofile_for_sflag,*outfile;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"h"))!=-1) {
    switch(c) {
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  argc-=optind;
  argv+=optind;

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  preofilename  = *argv;
  bd8filename = *++argv;
  coofilename = *++argv;
  coofilename_for_sflag = *++argv;
  outfilename = *++argv;

  read_ECEPE_parm(preofilename,bd8filename,&ECEPE_p,&nb_p);

  pi=acos(-1.0);

  bp=(int *)gcemalloc(sizeof(int)*(ECEPE_p.NUMATM-1)*2);
  bp_f=(int **)gcemalloc(sizeof(int *)*ECEPE_p.NUMATM);
  pair1_5=(int **)gcemalloc(sizeof(int *)*(ECEPE_p.NUMATM-1));
  pair1_4=(int **)gcemalloc(sizeof(int *)*(ECEPE_p.NUMATM-1));
  num1_4=(int *)gcemalloc(sizeof(int )*ECEPE_p.NUMATM);
  num1_5=(int *)gcemalloc(sizeof(int )*ECEPE_p.NUMATM);

  numb=(int *)gcemalloc(sizeof(int)*ECEPE_p.NUMATM);
  make_bd_pair_list(ECEPE_p,bp,bp_f,numb);

  make_dihed_pairs_list(ECEPE_p,bp_f,numb);

  make_int_pair_list(bp_f,numb,ECEPE_p.NUMATM,ECEPE_p.NUMATM,pair1_5,pair1_4,num1_4,num1_5);

  numnb=0;
  num14=0;
  for (i=0;i<ECEPE_p.NUMATM;++i)
    if ( num1_5[i]-1 > 0 )
      numnb+=num1_5[i]-1;
  for (i=0;i<ECEPE_p.NUMATM;++i)
    if ( num1_4[i]-1 > 0)
      num14+=num1_4[i]-1;

  pairs=(int *)gcemalloc(sizeof(int)*(numnb+num14)*2);
  k=0;
  for (i=0;i<ECEPE_p.NUMATM;++i) {
    for (j=0;j<num1_5[i]-1;++j) {
      pairs[k*2]=i+1;
      pairs[k*2+1]=pair1_5[i][j]+1;
      ++k;
    }
  }
  l=0;
  for (i=0;i<ECEPE_p.NUMATM;++i) {
    for(j=0;j<num1_4[i]-1;++j) {
      //    for(j=0;j<num1_4[i];++j) {
      pairs[(k+l)*2]=-(i+1);
      pairs[(k+l)*2+1]=pair1_4[i][j]+1;
      ++l;
    }
  }

  for (i=0;i<ECEPE_p.NUMATM;++i) {
    num1_5[i]-=1;
  }

  co=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMATM)*3);
  dihed=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMVAR));
  delta_dihed=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMVAR));
  for (i=0;i<ECEPE_p.NUMVAR;++i) delta_dihed[i]=0.0;

  co_dummy=(double *)gcemalloc(sizeof(double)*40*3);
  dihed_dummy=(double *)gcemalloc(sizeof(double)*10);
  ene=(double *)gcemalloc(sizeof(double)*6);

  coofile_for_sflag=efopen(coofilename_for_sflag,"r");
  read_ECEPE_detail_coo_cyc(coofile_for_sflag,co_dummy,dihed_dummy,ene);
  for (i=0;i<ECEPE_p.NUMVAR;++i) {
    dihed_dummy[i]=dihed_dummy[i]*pi/180.0;
    if (dihed_dummy[i]<-pi)
      dihed_dummy[i]+=2.0*pi;
    else if (dihed_dummy[i]>pi)
      dihed_dummy[i]-=2.0*pi;
  }
  k=0;
  for (i=0;i<ECEPE_p.NUMATM;++i) {
    for (j=0;j<3;++j) {
      co[i*3+j]=co_dummy[k];
      ++k;
    }
  }    
  calc_TORS_for_get_sabun(ECEPE_p.NUMVAR,co,ECEPE_p,dihed_dummy,delta_dihed);
  fclose(coofile_for_sflag);

  /**********************************************/
  /* for (i=0;i<ECEPE_p.NUMATM;++i)	        */
  /*   for (j=0;j<3;++j) 		        */
  /*     co[i*3+j]=ECEPE_p.atom[i].refcoord[j]; */
  /**********************************************/

  coofile=efopen(coofilename,"r");
  //  read_ECEPE_detail_coo(coofile,co,ECEPE_p.NUMATM);
  read_ECEPE_coo(coofile,co,dihed,ECEPE_p.NUMATM);
  fclose(coofile);

  outfile=efopen(outfilename,"w");
  calc_TORS_ECEPE2_for_check_out(ECEPE_p.NUMVAR,co,ECEPE_p,outfile,delta_dihed);
  fclose(outfile);

  return 0;
}

void USAGE(char *progname) {
  printf("-h -- help\n");
  printf("USAGE: %s profilename bd8filename coofilename coofilename_for_sflag outfilename\n", progname);
}

