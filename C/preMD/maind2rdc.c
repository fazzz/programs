#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <netcdf.h>
#include <ctype.h>
#include <getopt.h>

#include "ABAb.h"

#include "d2rc.h"

#include "PTL.h"
#include "EF.h"

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,m,o,ns,d;
  int numclut,nNumClutdummy,num,nNumClutLast;
  int flag,joinflag,flagsof,flagnum,flagINC;

  int nums,n;
  int numf;

  int typeflag=ATOM;
  int specifymove=OFF;
  int sidechainmove=ON;


  int numd2r;
  int numatom;
  int numres;
  int numresdummy;
  int **d2r;
  int *d2r_c;
  int da,db;
  int *resid;
  int *residdummy;

  CLTb *clt;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*clustfileinname,*parmfilename;
  char *clustfileoutname;

  FILE *inputfile,*clustfilein,*parmfile;
  FILE *clustfileout,*logfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"res",0,NULL,'e'},
    //    {"backbone",0,NULL,'b'},
    {"sidechain",0,NULL,'s'},
    {"move",0,NULL,'m'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hesm",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'e':
      typeflag=RESIDUE;
      break;
    case 'm':
      specifymove=ON;
      break;
    case 's':
      sidechainmove=OFF;
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
  inputfilename    = *argv;
  parmfilename  = *++argv;
  clustfileinname  = *++argv;
  clustfileoutname = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;

  clustfilein=efopen(clustfileinname,"r");
  clt=ABAbp_clustscan(clustfilein,&numclut);
  fclose(clustfilein);
  
  //  joinflag=0; for (i=0;i<numclut;++i) ABA_setJoin(clt,i,joinflag);
  clt[0].join=0;
  for(i=1; i<numclut; ++i) ABAb_setJoin(clt,i);

  inputfile=efopen(inputfilename,"r");
  resid=readd2rcinput(inputfile,&numres,specifymove);
  fclose(inputfile);

  if (sidechainmove==ON) {
    d2r=res2atom_BBc(resid,numatom,numres,&numd2r);
  }
  else {
    d2r=res2atom_SCc(resid,numres,&numd2r,clt,numclut);
  }

  for (i=0;i<numclut;++i) {
    clt[i].nNumClutOfChild=(int *)gcerealloc(clt[i].nNumClutOfChild,sizeof(int)*30);
    clt[i].terminal_atom_a=(int *)gcerealloc(clt[i].terminal_atom_a,sizeof(int)*30);
  }

  trans_d2rc(clt,&numclut,d2r,numd2r);

  clustfileout=efopen(clustfileoutname,"w");
  writed2rcoutput(clustfileout,clt,numclut);
  fclose(clustfileout);

  logfile=efopen("d2r_log.txt","w");
  for (i=0;i<numd2r;++i) {
    fprintf(logfile,"%5d - %5d\n",d2r[i][0],d2r[i][1]);
  }
  fclose(logfile);

  return 0;
  
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-m] specify move residue \n");
  printf("[-s] constraint side chain \n");
  printf("[-e] help \n");
  printf("%s [-h] inputfilename clustfileinname topfilename clustfileoutputfilename \n",progname);
}

