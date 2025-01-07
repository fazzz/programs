
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "FASYS.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int numstep;
  int flag=OFF;

  int sl=0;

  double v,q[5][3];

  char *line;
  size_t len=0;

  char *progname;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;


  double kd[2],n[2];

  char *inifilename,*outputfilename;
  FILE *inifile,*outputfile;

  kd[0]=2.5;
  kd[1]=2.5;
  n[0]=3.0;
  n[1]=3.0;

  progname=argv[0];

  while((c=getopt(argc,argv,"hfs:k:l:n:m:"))!=-1) {
    switch(c) {
    case 's':
      sl=atoi(optarg);
      break;
    case 'k':
      kd[0]=atof(optarg);
      break;
    case 'l':
      kd[1]=atof(optarg);
      break;
    case 'n':
      n[0]=atof(optarg);
      break;
    case 'm':
      n[1]=atof(optarg);
      break;
    case 'f':
      flag=ON;
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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  numstep = atoi(*argv);
  inifilename  = *++argv;
  outputfilename=*++argv;

  inifile=efopen(inifilename,"r");
  outputfile=efopen(outputfilename,"w"); 
  for (i=0;i<numstep;++i) {
    for (j=0;j<5;++j)
      for (k=0;k<3;++k)
	fscanf(inifile,"%lf",&q[j][k]);	
    v=FASYS_calcpote(q,kd,n);
    if (flag==ON)
      fprintf(outputfile,"%d %10.4lf\n",i,v);
    else
      fprintf(outputfile,"%10.4lf\n",v);
  }
  fclose(inifile);
  fclose(outputfile);

}

int USAGE(char *progname) {
  printf("USAGE:%s \n",progname);
  printf("-s sl (dissmiss lines)    \n");
  printf("-f formal output \n");
  printf("-k kd1 \n");
  printf("-l kd2 \n");
  printf("-m n1 \n");
  printf("-n n2 \n");
  printf("-h help  \n");
  printf("  numstep inifilename outputfilename\n");
}

 
