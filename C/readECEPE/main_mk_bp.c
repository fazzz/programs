#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ECEPE.h"
#include "PT.h"

#include "EF.h"


void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;

  struct ECEPE_parms pa;
  struct pnb nb_p;
  int *bp;
  int **bp_f;
  int **pair1_5,**pair1_4;
  int *num14,*num1_5;
  int *numb;

  char *progname;
  char *preofilename,*bd8filename;
  FILE *preofile,*bd8file;

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

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  preofilename  = *argv;
  bd8filename   = *++argv;

  read_ECEPE_parm(preofilename,bd8filename,&pa,&nb_p);

  bp=(int *)gcemalloc(sizeof(int)*(pa.NUMATM-1)*2);
  bp_f=(int **)gcemalloc(sizeof(int *)*pa.NUMATM);
  pair1_5=(int **)gcemalloc(sizeof(int *)*(pa.NUMATM-1));
  pair1_4=(int **)gcemalloc(sizeof(int *)*(pa.NUMATM-1));
  num14=(int *)gcemalloc(sizeof(int )*pa.NUMATM);
  num1_5=(int *)gcemalloc(sizeof(int )*pa.NUMATM);

  numb=(int *)gcemalloc(sizeof(int)*pa.NUMATM);
  make_bd_pair_list(pa,bp,bp_f,numb);

  make_int_pair_list(bp_f,numb,pa.NUMATM,pair1_5,pair1_4,num14,num1_5);

  /*******************************************************/
  /* for (i=0;i<pa.NUMATM;++i) {			 */
  /*   for (j=0;j<numb[i];++j) {			 */
  /*     printf("atom1=%d atom2=%d\n",i+1,bp_f[i][j]+1); */
  /*   }						 */
  /* }							 */
  /* printf("\n");					 */
  /*******************************************************/

  for (i=0;i<pa.NUMATM;++i) {
    printf("%d to ",i+1);
    for (j=0;j<num14[i];++j) {
      printf("%d ",pair1_4[i][j]+1);
    }
    printf("\n");
  }
  printf("\n");
  
  for (i=0;i<pa.NUMATM;++i) {
    printf(" %d to ",i+1);
    for (j=0;j<num1_5[i]-1;++j) {
      printf(" %d",pair1_5[i][j]+1);
    }
    printf("\n");
  }
  printf("\n");

  return 0;
}

void USAGE(char *progname) {
  printf("-h -- help\n");
  printf("USAGE: %s profilename bd8filename\n", progname);
}

