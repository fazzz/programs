#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ECEPE.h"
#include "EF.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k;
  int num;

  struct ECEPE_parms ECEPE_p;
  struct pnb nb_p;


  char *progname;

  char *preofilename;
  char *bd8filename;
  char *clustfilename;
  FILE *clustfile;

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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  preofilename  = *argv;
  bd8filename   = *++argv;
  clustfilename  = *++argv;

  read_ECEPE_parm(preofilename,bd8filename,&ECEPE_p,&nb_p);

  clustfile=efopen(clustfilename,"w");

  printf("%4d  \n",ECEPE_p.NUMVAR+1);
  fprintf(clustfile,"%4d  \n",ECEPE_p.NUMVAR+1);
  
  printf("%4d  ",ECEPE_p.atom[ECEPE_p.dihed[0].ibnd1-1].katom);
  fprintf(clustfile,"%4d  ",ECEPE_p.atom[ECEPE_p.dihed[0].ibnd1-1].katom);
  for (i=0;i<ECEPE_p.NUMVAR;++i) {
    printf("%4d  ",ECEPE_p.atom[ECEPE_p.dihed[i].ibnd2-1].katom);
    fprintf(clustfile,"%4d  ",ECEPE_p.atom[ECEPE_p.dihed[i].ibnd2-1].katom);
  }
  printf("  \n");
  fprintf(clustfile,"  \n");

  printf("   1 ");
  fprintf(clustfile,"   1 ");
  for (i=0;i<ECEPE_p.NUMVAR;++i) {
    printf("   1 ");
    fprintf(clustfile,"   1 ");
  }
  printf("  \n");
  fprintf(clustfile,"  \n");

  printf("   4   2 ");
  fprintf(clustfile,"   4   2 ");
  for (i=0;i<ECEPE_p.NUMVAR+1-4;++i) {
    printf(" ");
  }
  printf("   2   4 ");
  fprintf(clustfile,"   2    4 ");
  printf("  \n");
  fprintf(clustfile,"  \n");

  printf("   1 ");
  fprintf(clustfile,"   1 ");
  for (i=0;i<ECEPE_p.NUMVAR;++i) {
    printf("   1 ");
    fprintf(clustfile,"   1 ");
  }
  printf("  \n");
  fprintf(clustfile,"  \n");

  for (i=0;i<ECEPE_p.NUMVAR+1;++i) {
    printf("   3 ");
    fprintf(clustfile,"   3 ");
  }
  printf("  \n");
  fprintf(clustfile,"  \n");

  printf("%4d  ",ECEPE_p.atom[ECEPE_p.dihed[0].ibnd1-1].katom);
  fprintf(clustfile,"%4d  ",ECEPE_p.atom[ECEPE_p.dihed[0].ibnd1-1].katom);
  for (i=0;i<ECEPE_p.NUMVAR;++i) {
    printf("%4d  ",ECEPE_p.atom[ECEPE_p.dihed[i].ibnd2-1].katom);
    fprintf(clustfile,"%4d  ",ECEPE_p.atom[ECEPE_p.dihed[i].ibnd2-1].katom);
  }
  printf("  \n");
  fprintf(clustfile,"  \n");

  printf("   0 ");
  fprintf(clustfile,"   0 ");
  for (i=0;i<ECEPE_p.NUMVAR;++i) {
    ;
  }
  printf("  \n");
  fprintf(clustfile,"  \n");

  printf("   2 ");
  fprintf(clustfile,"   2 ");
  for (i=0;i<ECEPE_p.NUMVAR+1-2;++i) {    
    printf("     ");
    fprintf(clustfile,"     ");
  }
  printf("  -1 ");
  fprintf(clustfile,"  -1 ");
  printf("  \n");
  fprintf(clustfile,"  \n");

  for (i=0;i<ECEPE_p.NUMVAR+1;++i) {
    printf("%4d ",i+1);
    fprintf(clustfile,"%4d ",i+1);
  }
  printf("  \n");
  fprintf(clustfile,"  \n");

  printf("0 \n");
  fprintf(clustfile,"0 \n");

  fclose(clustfile);
  
  return 0;
}

void USAGE(char *progname) {
  printf("-h -- help\n");
  printf("USAGE: %s profilename bd8filename massfilename\n", progname);
}

