#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <netcdf.h>
#include <ctype.h>
#include <getopt.h>

#include "ABA_hosoku.h"
//#include "ABA.h"

#include "d2r.h"

#include "PTL.h"
#include "EF.h"

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,m,o,ns,d;
  int numclut,numatom;

  int num=100;

  CLTh *clt;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *clustfileinname,*parmfilename;
  char *clustfileoutname;

  FILE *clustfilein,*parmfile;
  FILE *clustfileout,*logfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"num",1,NULL,'n'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hn:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'n':
      num=atoi(optarg);
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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  parmfilename     = *argv;
  clustfileinname  = *++argv;
  clustfileoutname = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;

  clustfilein=efopen(clustfileinname,"r");
  clt=ABAhp_clustscan(clustfilein,&numclut);
  fclose(clustfilein);

  numclut=num;

  clt[numclut-1].terminal=0;
  clt[numclut-1].num_branch=1;
  clt[numclut-1].nNumClutOfChild[0]=-1;
  
  clustfileout=efopen(clustfileoutname,"w");
  writed2routput(clustfileout,clt,numclut);
  fclose(clustfileout);

  return 0;
  
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[--num] number of clust \n");
  printf("%s [-h] clustfileinname topfilename clustfileoutputfilename \n",progname);
}

