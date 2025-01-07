#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <getopt.h>

#include "PTL.h"
#include "EF.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,m;
  int nl[20];
  int dtype;
  int d[4];
  double p[20],f[20];
  int numatom,numres;
  double pi;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*inputfilename2,*parmfilename,*name;
  char names[100][6],temp[5];

  FILE *inputfile,*inputfile2,*parmfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"h",long_opt,&opt_idx))!=-1) {
    switch(c) {
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
  inputfilename     = *argv;
  inputfilename2    = *++argv;
  parmfilename      = *++argv;
  name = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  
  inputfile=efopen(inputfilename,"r");
  for (i=0;i<20;++i) {
    fscanf(inputfile,"%6s",&names[i]);
    fscanf(inputfile,"%d",&nl[i]);
    fscanf(inputfile,"%lf",&p[i]);
    fscanf(inputfile,"%lf",&f[i]);
  }
  fclose(inputfile);

  for (i=0;i<20;++i) {
    if (strncmp(name,names[i],3)==0)
      break;
  }
  if (i==20) {
    printf("error\n");
    exit(1);
  }
  m=i;

  inputfile2=efopen(inputfilename2,"r");
  for (i=0;i<nl[m]-1;++i)  getline(&line,&len,inputfile2);
  fscanf(inputfile2,"%5s",&temp);
  fscanf(inputfile2,"%d",&d[0]);
  fscanf(inputfile2,"%d",&d[1]);
  fscanf(inputfile2,"%d",&d[2]);
  fscanf(inputfile2,"%d",&d[3]);
  fclose(inputfile2);

  for (i=0;i<AP.NPHIH;++i) {
    if ((abs(AP.PH[i][1])/3+1)==d[1] && (abs(AP.PH[i][2])/3+1)==d[2]) {
      dtype = AP.PH[i][4]-1;
      
      printf("%3d-%3d-%3d-%3d %4s-%4s-%4s-%4s F=%3e N=%3d PHASE=%3lf type=%3d\n",
	     (abs(AP.PH[i][0]))/3+1,(abs(AP.PH[i][1]))/3+1,(abs(AP.PH[i][2]))/3+1,(abs(AP.PH[i][3]))/3+1,
	     AP.IGRAPH[(abs(AP.PH[i][0]))/3],AP.IGRAPH[(abs(AP.PH[i][1]))/3],
	     AP.IGRAPH[(abs(AP.PH[i][2]))/3],AP.IGRAPH[(abs(AP.PH[i][3]))/3],
	     AP.PK[dtype],AP.PN[dtype],AP.PHASE[dtype],dtype+1);
    }
  }

  for (i=0;i<AP.MPHIA;++i) {
    if ((abs(AP.PA[i][1])/3+1)==d[1] && (abs(AP.PA[i][2])/3+1)==d[2]) {
      dtype = AP.PA[i][4]-1;
      
      printf("%3d-%3d-%3d-%3d %4s-%4s-%4s-%4s F=%3e N=%3d PHASE=%3lf type=%3d\n",
	     (abs(AP.PA[i][0]))/3+1,(abs(AP.PA[i][1]))/3+1,(abs(AP.PA[i][2]))/3+1,(abs(AP.PA[i][3]))/3+1,
	     AP.IGRAPH[(abs(AP.PA[i][0]))/3],AP.IGRAPH[(abs(AP.PA[i][1]))/3],
	     AP.IGRAPH[(abs(AP.PA[i][2]))/3],AP.IGRAPH[(abs(AP.PA[i][3]))/3],
	     AP.PK[dtype],AP.PN[dtype],AP.PHASE[dtype],dtype+1);
    }
  }

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename refcrdfilename parmfilename outputfilename outputfilename2 trjfilename\n",progname);
}
