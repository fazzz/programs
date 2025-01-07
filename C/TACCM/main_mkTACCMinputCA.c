#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PTL.h"
#include "FFL.h"
#include "EF.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,num;
  int numatom;

  int CAflag=OFF;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *TACCMfilename,*parmfilename;
  FILE *TACCMfile,*parmfile;

  char *progname;
  int opt_idx=1;

  struct option long_opt[] = {
    {"h",0,NULL,'h'},
    {"CA",0,NULL,'a'},
    {"HA",0,NULL,'H'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"haH",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'a':
      CAflag=ON;
      break;
    case 'H':
      CAflag=OFF;
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

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  parmfilename   = *argv;
  TACCMfilename  = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;

  num=0;
  for (i=0;i<numatom;++i) {
    if (CAflag==ON) {
      if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
	++num;
      }
    }
    else {
      if (strncmp(AP.IGRAPH[i],"H",1)!=0) {
	++num;
      }
    }
  }

  TACCMfile=efopen(TACCMfilename,"w");
  fprintf(TACCMfile,"%4d\n",num);
  for (i=0;i<numatom;++i) {
    if (CAflag==ON) {
      if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
	fprintf(TACCMfile,"%4d \n",i+1);
      }
    }
    else {
      if (strncmp(AP.IGRAPH[i],"H",1)!=0) {
	fprintf(TACCMfile,"%4d \n",i+1);
      }
    }
  }
  fclose(TACCMfile);

  return 0;

}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-CA] CAflag \n");
  printf("[-HA] HAflag \n");
  printf("[-h] help \n");
  printf("%s [-h] parmfilename TACCMfilename \n",progname);
}
