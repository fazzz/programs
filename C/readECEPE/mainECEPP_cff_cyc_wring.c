#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ECEPE.h"
#include "PROTOPO.h"
#include "PT.h"
#include "FF.h"

#include "EF.h"

#define ON  1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int d;
  int checkflag=OFF,sflag=OFF,fflag=OFF,aromflag=ON;

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

  int **matrix;

  char *name_atom_list;

  struct ECEPE_parms ECEPE_p;
  struct pnb nb_p;
  struct ECEPE_pote p;
  struct ECEPE_force f;

  char *progname;
  char *preofilename,*bd8filename,*parmtopfilename,*indexfilename,*coofilename,*dfilename,*efilename,*coofilename_for_sflag;
  FILE *parmtopfile,*indexfile,*coofile,*dp,*dfile,*efile,*coofile_for_sflag;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"has"))!=-1) {
    switch(c) {
    case 'a':
      checkflag=ON;
      break;
    case 's':
      sflag=ON;
      break;
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

  if (argc < /*2*/5) {
    USAGE(progname);
    exit(1);
  }
  preofilename  = *argv;
  bd8filename = *++argv;
  coofilename = *++argv;
  dfilename = *++argv;
  efilename = *++argv;
  if (sflag == ON)
    coofilename_for_sflag = *++argv;

  read_ECEPE_parm(preofilename,bd8filename,&ECEPE_p,&nb_p);

  pi=acos(-1.0);

  bp=(int *)gcemalloc(sizeof(int)*(ECEPE_p.NUMATM-1)*2);
  bp_f=(int **)gcemalloc(sizeof(int *)*ECEPE_p.NUMATM);
  pair1_5=(int **)gcemalloc(sizeof(int *)*(ECEPE_p.NUMATM-1));
  pair1_4=(int **)gcemalloc(sizeof(int *)*(ECEPE_p.NUMATM-1));
  num1_4=(int *)gcemalloc(sizeof(int )*ECEPE_p.NUMATM);
  num1_5=(int *)gcemalloc(sizeof(int )*ECEPE_p.NUMATM);

  numb=(int *)gcemalloc(sizeof(int)*ECEPE_p.NUMATM);
  //  make_bd_pair_list(ECEPE_p,bp,bp_f,numb);
  name_atom_list=(char *)gcemalloc(sizeof(char)*ECEPE_p.NUMATM*4);
  for (i=0;i<ECEPE_p.NUMATM;++i) for (j=0;j<4;++j)
	name_atom_list[i*4+j]=ECEPE_p.atom[i].name_atom[j];
  make_bp(ECEPE_p.NUMATM,bp,bp_f,numb,name_atom_list,aromflag);

  make_dihed_pairs_list(ECEPE_p,bp_f,numb);

  matrix=(int **)gcemalloc(sizeof(int *)*ECEPE_p.NUMATM);
  for (i=0;i<ECEPE_p.NUMATM;++i) matrix[i]=(int *)gcemalloc(sizeof(int)*ECEPE_p.NUMATM);

  make_nb_matrix(bp_f,numb,3,matrix,ECEPE_p.NUMATM);

  set_nb_pairs(matrix,ECEPE_p.NUMATM,pair1_5,pair1_4,num1_5,num1_4);

  //  make_int_pair_list(bp_f,numb,ECEPE_p.NUMATM,ECEPE_p.NUMATM,pair1_5,pair1_4,num1_4,num1_5);

  numnb=0;
  num14=0;
  for (i=0;i<ECEPE_p.NUMATM;++i) numnb+=num1_5[i];
  for (i=0;i<ECEPE_p.NUMATM;++i) num14+=num1_4[i];

  pairs=(int *)gcemalloc(sizeof(int)*(numnb+num14)*2);
  k=0;
  for (i=0;i<ECEPE_p.NUMATM;++i) {
    for (j=0;j<num1_5[i];++j) {
      pairs[k*2]=i+1;
      pairs[k*2+1]=pair1_5[i][j]+1;
      ++k;
    }
  }
  l=0;
  for (i=0;i<ECEPE_p.NUMATM;++i) {
    for(j=0;j<num1_4[i];++j) {
      pairs[(k+l)*2]=-(i+1);
      pairs[(k+l)*2+1]=pair1_4[i][j]+1;
      ++l;
    }
  }

  /************************************/
  /* for (i=0;i<ECEPE_p.NUMATM;++i) { */
  /*   num1_5[i]-=1;		      */
  /* }				      */
  /************************************/

  co=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMATM)*3);
  dihed=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMVAR));
  delta_dihed=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMVAR));
  for (i=0;i<ECEPE_p.NUMVAR;++i) delta_dihed[i]=0.0;

  co_dummy=(double *)gcemalloc(sizeof(double)*40*3);
  dihed_dummy=(double *)gcemalloc(sizeof(double)*10);
  ene=(double *)gcemalloc(sizeof(double)*6);

  if ( sflag == ON ) {
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
  }

  coofile=efopen(coofilename,"r");
  dfile=efopen(dfilename,"w");
  efile=efopen(efilename,"w");

  while ((c=read_ECEPE_detail_coo_cyc(coofile,co_dummy,dihed_dummy,ene))!=-1) {
    k=0;
    for (i=0;i<ECEPE_p.NUMATM;++i) {
      for (j=0;j<3;++j) {
	ECEPE_p.atom[i].refcoord[j]=co_dummy[k];
	++k;
      }
    }    

    if (checkflag==ON) {
      for (i=0;i<4;++i) {
	printf("atom num %d=",i);
	scanf("%d",&atom[i]);
	printf("\n");
      }
    
      dihed_check=calc_TORS_for_check(co,atom[0],atom[1],atom[2],atom[3]);
    
      printf("%lf \n",dihed_check*180/pi);
    }

    /*************************************************************************/
    /* printf("-----------------------------------\n");			     */
    /* for (i=0;i<ECEPE_p.NUMATM;++i) {					     */
    /*   for (j=0;j<4;++j) printf("%c",ECEPE_p.atom[i].name_atom[j]);	     */
    /*   printf("%d -- ",num1_5[i]);					     */
    /*   for (j=0;j<num1_5[i];++j) {					     */
    /* 	for (k=0;k<4;++k)						     */
    /* 	  printf("%c",ECEPE_p.atom[(pair1_5[i][j])].name_atom[k]);	     */
    /* 	printf("(%d)",pair1_5[i][j]+1);					     */
    /*   }								     */
    /*   printf("\n");							     */
    /* }								     */
    /* printf("\n");							     */
    /* for (i=0;i<ECEPE_p.NUMATM;++i) {					     */
    /*   for (j=0;j<4;++j) printf("%c",ECEPE_p.atom[i].name_atom[j]);	     */
    /*   printf("%d -- ",num1_4[i]);					     */
    /*   for (j=0;j<num1_4[i];++j) {					     */
    /* 	for (k=0;k<4;++k)						     */
    /* 	  printf("%c",ECEPE_p.atom[(pair1_4[i][j])].name_atom[k]);	     */
    /* 	printf("(%d)",pair1_4[i][j]+1);					     */
    /*   }								     */
    /*   printf("\n");							     */
    /* }								     */
    /* printf("-----------------------------------\n");			     */
    /*************************************************************************/

    calc_ff_ECEPE_for_db_cyc(&p,&f,ECEPE_p,nb_p,pair1_5,pair1_4,num1_5,num1_4,dfile,delta_dihed);
  
    fprintf(efile,"%14.10e ",p.p_tors);
    fprintf(efile,"%14.10e ",p.p_es);
    fprintf(efile,"%14.10e ",p.p_nb);
    fprintf(efile,"   0.0000000000e+00 ");
    fprintf(efile,"   0.0000000000e+00 ");
    fprintf(efile,"%14.10e \n",p.p_t);
  }

  fclose(coofile);
  fclose(dfile);
  fclose(efile);

  return 0;
}

void USAGE(char *progname) {
  printf("-a -- check flag for dihedral anle\n");
  printf("-s -- sabun to dihed in inspidas\n");
  printf("-h -- help\n");
  printf("USAGE: %s profilename bd8filename coofilename dfilename efilename\n", progname);
}

