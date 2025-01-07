#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ECEPE.h"
#include "EF.h"

#define ON  1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j;
  int d1,d2;
  int dt1,dt2;
  int na,nb;

  int numatom;

  int numtotal1_4,numtotal1_5;
  int **pairs1_5,**pairs1_4;
  int *num1_4,*num1_5;
  int **matrixpairs;

  struct ECEPE_parms ECEPE_p;
  struct pnb nb_p;
  struct ECEPE_pote p;
  struct ECEPE_force f;

  char *progname;
  char *preofilename,*bd8filename;
  char *inputpairsfilename,*outputpairsfilename;
  FILE *inputpairsfile,*outputpairsfile;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"has"))!=-1) {
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  preofilename  = *argv;
  bd8filename = *++argv;
  inputpairsfilename = *++argv;
  outputpairsfilename = *++argv;

  read_ECEPE_parm(preofilename,bd8filename,&ECEPE_p,&nb_p);
  numatom=ECEPE_p.NUMATM;
  pairs1_4=(int **)gcemalloc(sizeof(int *)*numatom);
  pairs1_5=(int **)gcemalloc(sizeof(int *)*numatom);
  num1_4=(int *)gcemalloc(sizeof(int)*numatom);
  num1_5=(int *)gcemalloc(sizeof(int)*numatom);
  matrixpairs=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) matrixpairs[i]=(int *)gcemalloc(sizeof(int)*numatom);

  for (i=0;i<numatom;++i) {
    num1_4[i]=0;
    num1_5[i]=0;
    for (j=0;j<numatom;++j) matrixpairs[i][j]=-1;
  }

  inputpairsfile=efopen(inputpairsfilename,"r");
  fscanf(inputpairsfile,"%d",&numtotal1_4);
  fscanf(inputpairsfile,"%d",&numtotal1_5);
  for (i=0;i<numtotal1_4;++i) {
    fscanf(inputpairsfile,"%d",&d1);
    fscanf(inputpairsfile,"%d",&d2);
    dt1=ECEPE_p.atom[d1-1].katom;
    dt2=ECEPE_p.atom[d2-1].katom;
    matrixpairs[dt1-1][dt2-1]=3;
    matrixpairs[dt2-1][dt1-1]=3;
  }
  for (i=0;i<numtotal1_5;++i) {
    fscanf(inputpairsfile,"%d",&d1);
    fscanf(inputpairsfile,"%d",&d2); 
    dt1=ECEPE_p.atom[d1-1].katom;
    dt2=ECEPE_p.atom[d2-1].katom;
    matrixpairs[dt1-1][dt2-1]=4;
    matrixpairs[dt2-1][dt1-1]=4;
  }
  fclose(inputpairsfile);

  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (matrixpairs[i][j]==3) {
	pairs1_4[i]=(int *)gcerealloc(pairs1_4[i],sizeof(int)*(num1_4[i]+1));
	pairs1_4[i][num1_4[i]]=j;
	num1_4[i]+=1;
      }
      else if (matrixpairs[i][j]==4) {
	pairs1_5[i]=(int *)gcerealloc(pairs1_5[i],sizeof(int)*(num1_5[i]+1));
	pairs1_5[i][num1_5[i]]=j;
	num1_5[i]+=1;
      }
    }
  }

  numtotal1_4=0;
  for (i=0;i<numatom;++i) numtotal1_4+=num1_4[i];
  numtotal1_5=0;
  for (i=0;i<numatom;++i) numtotal1_5+=num1_5[i];

  outputpairsfile=efopen(outputpairsfilename,"w");
  fprintf(outputpairsfile,"%d\n",numtotal1_4);
  fprintf(outputpairsfile,"%d\n",numtotal1_5);
  for (i=0;i<numatom;++i) {
    for (j=0;j<num1_4[i];++j) {
      fprintf(outputpairsfile,"%d ",i+1);
      fprintf(outputpairsfile,"%d\n",pairs1_4[i][j]+1);
    }
  }
  for (i=0;i<numatom;++i) {
    for (j=0;j<num1_5[i];++j) {
      fprintf(outputpairsfile,"%d ",i+1);
      fprintf(outputpairsfile,"%d\n",pairs1_5[i][j]+1);
    }
  }
  fclose(outputpairsfile);

}

void USAGE(char *progname) {
  printf("-h -- help\n");
  printf("USAGE: %s profilename bd8filename inputpairsfilename outputpairsfilename\n", progname);
}
