#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ECEPE.h"
#include "EF.h"

#define ARRAY_SIZE(array) (sizeof(array)/sizeof(array[0]))

#define ON  1
#define OFF 0

struct nb_interact {
 int i;
 int j;
   double ene;
};

int compare_nb(struct nb_interact*a , struct nb_interact*b);

void USAGE(char *progname);


int main(int argc, char *argv[]) {
  int i,j,n;
  int d1,d2;
  double d3;
  int numatom;
  int numline;

  int *indexjtok;
  int *indexktoj;

  struct nb_interact *nb_int;

  struct ECEPE_parms ECEPE_p;
  struct pnb nb_p;

  char *progname;
  char *preofilename,*bd8filename;
  char *inputfilename,*outputfilename;
  FILE *inputfile,*outputfile;

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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  numline = atoi(*argv);
  preofilename  = *++argv;
  bd8filename = *++argv;
  inputfilename = *++argv;
  outputfilename = *++argv;

  read_ECEPE_parm(preofilename,bd8filename,&ECEPE_p,&nb_p);
  numatom=ECEPE_p.NUMATM;
  indexjtok=(int *)gcemalloc(sizeof(int)*numatom);
  indexktoj=(int *)gcemalloc(sizeof(int)*numatom);
  for (i=0;i<numatom;++i) {
    indexjtok[ECEPE_p.atom[i].jatom-1]=ECEPE_p.atom[i].katom;
    indexktoj[ECEPE_p.atom[i].katom-1]=ECEPE_p.atom[i].jatom;
  }

  nb_int=(struct nb_interact *)gcemalloc(sizeof(struct nb_interact)*numline);
  inputfile=efopen(inputfilename,"r");
  for (i=0;(c=fscanf(inputfile,"%d",&d1))!=-1;++i) {
    fscanf(inputfile,"%d",&d2);
    fscanf(inputfile,"%lf",&d3);
    //    nb_int=(int *)gcerealloc(nb_int,sizeof(int)*(i+1)*2);
    nb_int[i].i=indexjtok[d1]-1;
    nb_int[i].j=indexjtok[d2]-1;
    nb_int[i].ene=d3;
  }
  fclose(inputfile);
  n=i;

  qsort(nb_int,ARRAY_SIZE(nb_int),sizeof(struct nb_interact),compare_nb);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<n;++i) {
    fprintf(outputfile,"%d %d %lf\n",nb_int[i].i+1,nb_int[i].j+1,nb_int[i].ene);
  }
  fclose(outputfile);

}

void USAGE(char *progname) {
   printf("-h -- help\n");
   printf("USAGE: %s profilename bd8filename inputfilename outputfilename\n", progname);
 }

int compare_nb(struct nb_interact*a , struct nb_interact*b) {
  int ia,ib;
  int ja,jb;

  ia=a->i;
  ib=b->i;

  if ((int*)ia==(int*)ib) {
    ja=a->j;
    jb=b->j;
    return (*(int*)ja -*(int*)jb);
  }
  else{
    return (*(int*)ia -*(int*)ib);
  }
}
