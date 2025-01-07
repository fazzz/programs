
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "FF.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int flag;

  char *line;
  size_t len=0;

  char *progname;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  int numindexnb,numindex14;
  int *indexnb,*index14;

  char *parmtopname,*logfilename;
  FILE *parmtop,*logfile;

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
  parmtopname  = *argv;
  logfilename  = *++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtop(parmtop);
  fclose(parmtop);

  ff_set_non_bonding_index_1(&numindexnb,&numindex14);
  indexnb=(int *)gcemalloc(sizeof(int)*numindexnb*2);
  index14=(int *)gcemalloc(sizeof(int)*numindex14*2);
  ff_set_non_bonding_index_2(indexnb,index14);

  logfile=efopen(logfilename,"w");
  for(i=0;i<numindexnb;++i) fprintf(logfile,"%d - %d \n",indexnb[i*2]+1,indexnb[i*2+1]+1);
  fprintf(logfile,"\n");
  for(i=0;i<numindex14;++i) fprintf(logfile,"%d - %d \n",index14[i*2]+1,index14[i*2+1]+1);
  fprintf(logfile,"\n  ");
  fclose(logfile);

}

int USAGE(char *progname) {
  printf("USAGE:%s \n",progname);
  printf("-t dt    \n");
  printf("-h help  \n");
  printf("  parmtopname\n");
}

 
