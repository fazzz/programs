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
  int cofileflag=OFF,checkflag=OFF;

  int numnb,num14;
  int *indexnb,*index14,*pairs;

  double *co;
  double *dihed;
  double pi;
  double dihed_check;

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
  char *preofilename,*bd8filename,*parmtopfilename,*indexfilename,*coofilename;
  FILE *parmtopfile,*indexfile,*coofile,*dp;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"hca"))!=-1) {
    switch(c) {
    case 'a':
      checkflag=ON;
      break;
    case 'c':
      cofileflag=ON;
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

  if (argc < /*2*/3) {
    USAGE(progname);
    exit(1);
  }
  preofilename  = *argv;
  bd8filename = *++argv;
  //  parmtopfilename = *++argv;
  //  indexfilename = *++argv;
  coofilename = *++argv;

  read_ECEPE_parm(preofilename,bd8filename,&ECEPE_p,&nb_p);
  /********************************************/
  /* parmtopfile=efopen(parmtopfilename,"r"); */
  /* readParmtop(parmtopfile);		      */
  /* fclose(parmtopfile);		      */
  /********************************************/

  /****************************************/
  /* printf("NATOM=%4d ",ECEPE_p.NUMATM); */
  /* printf("NDIHE=%4d ",ECEPE_p.NUMVAR); */
  /* printf("NRESD=%4d ",ECEPE_p.NUMRES); */
  /* printf("NINTE=%4d ",ECEPE_p.NUMINT); */
  /* printf("NSS  =4%d\n",ECEPE_p.NUMS);  */
  /****************************************/

  /**********************************************/
  /* for (i=0;i<ECEPE_p.NUMVAR;++i) {	        */
  /*   printf("A=%10.4lf ",ECEPE_p.dihed[i].A); */
  /*   printf("NB=%4d ",ECEPE_p.dihed[i].NB);   */
  /*   printf("NS=%4d\n",ECEPE_p.dihed[i].NS);  */
  /* }					        */
  /**********************************************/

  /*******************************************************/
  /* for (i=0;i<ECEPE_p.NUMATM;++i) {			 */
  /*   for (j=0;j<3;++j)				 */
  /*     printf("%10.6lf ",ECEPE_p.atom[i].refcoord[j]); */
  /*   printf("charge=%10.6lf ",ECEPE_p.atom[i].charge); */
  /*   printf("nbtype=%4d \n",ECEPE_p.atom[i].nbtype);	 */
  /* }							 */
  /*******************************************************/

  /********************************************************/
  /* ff_set_non_bonding_index_1(&numnb,&num14);		  */
  /* indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);	  */
  /* index14=(int *)gcemalloc(sizeof(int)*num14*2);	  */
  /* ff_set_non_bonding_index_2(indexnb,index14);	  */
  /* pairs=(int *)gcemalloc(sizeof(int)*(numnb+num14)*2); */
  /* 							  */
  /* for (i=0;i<numnb;++i) {				  */
  /*   pairs[i*2]=indexnb[i*2];				  */
  /*   pairs[i*2+1]=indexnb[i*2+1];			  */
  /* }							  */
  /* for (i=numnb;i<numnb+num14;++i) {			  */
  /*   pairs[i*2]=-1.0*index14[(i-numnb)*2];		  */
  /*   pairs[i*2+1]=index14[(i-numnb)*2+1];		  */
  /* }							  */
  /********************************************************/

  /****************************************/
  /* j=0;				  */
  /* indexfile=efopen(indexfilename,"r"); */
  /* for (i=0;i<ECEPE_p.NUMATM;) {	  */
  /*   fscanf(indexfile,"%d",&d);	  */
  /*   if (d==0) ++i;			  */
  /*   else {				  */
  /*     pairs[j*2]=i+1;		  */
  /*     pairs[j*2+1]=d;		  */
  /*     ++j;				  */
  /*   }  				  */
  /* }					  */
  /* k=0;				  */
  /* for (i=0;i<ECEPE_p.NUMATM;) {	  */
  /*   fscanf(indexfile,"%d",&d);	  */
  /*   if (d==0) ++i;			  */
  /*   else {				  */
  /*     pairs[(k+j)*2]=-(i+1);		  */
  /*     pairs[(k+j)*2+1]=d;		  */
  /*     ++k;				  */
  /*   }  				  */
  /* }					  */
  /* fclose(indexfile);			  */
  /****************************************/

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

  /****************************************************/
  /* dp=efopen("dp","r");			      */
  /* for (i=0;i<ECEPE_p.NUMVAR;++i) {		      */
  /*   for (j=0;j<4;++j) {			      */
  /*     fscanf(dp,"%d",&ECEPE_p.dihed[i].dpairs[j]); */
  /*   }					      */
  /* }						      */
  /* fclose(dp);				      */
  /****************************************************/

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
    for(j=0;j<num1_4[i];++j) {
      pairs[(k+l)*2]=-(i+1);
      pairs[(k+l)*2+1]=pair1_4[i][j]+1;
      ++l;
    }
  }

  for (i=0;i<ECEPE_p.NUMATM;++i) {
    num1_5[i]-=1;
  }
  /**************************************/
  /* for (i=0;i<ECEPE_p.NUMATM;++i) {   */
  /*   printf("%d :", i+1);	        */
  /*   for (j=0;j<num1_5[i];++j) {      */
  /*     printf("%d ",pair1_5[i][j]+1); */
  /*   }			        */
  /*   printf("\n");		        */
  /* }				        */
  /**************************************/
  
  /**************************************/
  /* for (i=0;i<ECEPE_p.NUMATM;++i) {   */
  /*   printf("%d :", i+1);	        */
  /*   for (j=0;j<num1_4[i];++j) {      */
  /*     printf("%d ",pair1_4[i][j]+1); */
  /*   }			        */
  /*   printf("\n");		        */
  /* }				        */
  /**************************************/

  co=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMATM)*3);
  dihed=(double *)gcemalloc(sizeof(double)*(ECEPE_p.NUMVAR));

  coofile=efopen(coofilename,"r");
  if ( cofileflag == OFF )
  read_ECEPE_coo(coofile,co,dihed,ECEPE_p.NUMATM);
  else
    read_ECEPE_detail_coo(coofile,co,ECEPE_p.NUMATM);
  fclose(coofile);

  for (i=0;i<ECEPE_p.NUMATM;++i)
    for (j=0;j<3;++j)
      ECEPE_p.atom[i].refcoord[j]=co[i*3+j];

  for (i=0;i<ECEPE_p.NUMVAR;++i)
    ECEPE_p.dihed[i].angle=dihed[i]*pi/180.0;

  if (checkflag==ON) {
    for (i=0;i<4;++i) {
      printf("atom num %d=",i);
      scanf("%d",&atom[i]);
      printf("\n");
    }
    
    dihed_check=calc_TORS_for_check(co,atom[0],atom[1],atom[2],atom[3]);
    
    printf("%lf \n",dihed_check*180/pi);
  }

  calc_ff_ECEPE_for_db(&p,&f,ECEPE_p,nb_p,pair1_5,pair1_4,num1_5,num1_4/*k+j*//*(numnb+num14)*/);
  //  calc_ff_ECEPE(&p,&f,ECEPE_p,nb_p,pairs,/*k+j*/(numnb+num14));

  /****************************************************/
  /* for (i=0;i<k+j;++i) {			      */
  /*   printf("a=%d b=%d\n",pairs[i*2],pairs[i*2+1]); */
  /* }						      */
  /* printf("\n");				      */
  /****************************************************/
  
  printf("%10.6e ",p.p_tors);
  printf("%10.6e ",p.p_es);
  printf("%10.6e ",p.p_nb);
  printf("0.000000e+00 ");
  printf("0.000000e+00 ");
  printf("%10.6e ",p.p_t);

  return 0;
}

void USAGE(char *progname) {
  printf("-a -- check flag for dihedral anle\n");
  printf("-c -- detail coordinate file\n");
  printf("-h -- help\n");
  printf("USAGE: %s profilename bd8filename coofilename\n", progname);
}

