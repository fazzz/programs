
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "MB.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,sl=0;
  int flagKOP='K';

  double pi;
  
  char *line;
  size_t len=0;

  char *progname;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  int numatom,numstep,numsdihed,numdihedtotal;
  int *sdihednum;
  double *dihed,*sdihed,*dihed_seled;
  double *crd,*crd_seled;
  double mind=0.0,temp;
  
  char *trjfilename,*outputfilename,*parmtopfilename;
  FILE *trjfile,*outputfile,*parmfile;
  
  progname=argv[0];

  pi=acos(-1.0);
  
  while((c=getopt(argc,argv,"hkops:"))!=-1) {
    switch(c) {
    case 'p':
      flagKOP='P';
      break;
    case 'o':
      flagKOP='O';
      break;
    case 'k':
      flagKOP='K';
      break;
    case 's':
      sl=atoi(optarg);
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
  
  if (argc < 4 ) {
    USAGE(progname);
    exit(1);
  }
  numstep=atoi(*argv);
  numsdihed=atoi(*++argv);
  sdihednum=(int *)gcemalloc(sizeof(int)*numsdihed);
  sdihed=(double *)gcemalloc(sizeof(double)*numsdihed);
  dihed_seled=(double *)gcemalloc(sizeof(double)*numsdihed);
  for (i=0;i<numsdihed;++i) {
    sdihednum[i]=atoi(*++argv)-1;
    sdihed[i]=atof(*++argv)*pi/180.0;
  }
  trjfilename=*++argv;
  parmtopfilename=*++argv;
  outputfilename=*++argv;

  parmfile=efopen(parmtopfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  crd_seled=(double *)gcemalloc(sizeof(double)*numatom*3);

  numdihedtotal=countdihed(flagKOP);
  dihed=(double *)gcemalloc(sizeof(double)*numdihedtotal);
  
  trjfile=efopen(trjfilename,"r");
  for (i=0;i<sl;++i) getline(&line,&len,trjfile);
  for (i=0;i<numstep;++i) {
    io_scantrj2(trjfile,numatom,1,crd);
    dihed=CD(crd,flagKOP);
    temp=0.0;
    for (j=0;j<numsdihed;++j) temp+=(dihed[sdihednum[j]]-sdihed[j])*(dihed[sdihednum[j]]-sdihed[j]);
    temp=sqrt(temp);
    if (i>0) {
      if (mind>temp) {
	mind=temp;
	for (j=0;j<numatom*3;++j) crd_seled[j]=crd[j];
      }
    }
    else {
      mind=temp;
      for (j=0;j<numatom*3;++j) crd_seled[j]=crd[j];
    }
  }
  fclose(trjfile);
  
  outputfile=efopen(outputfilename,"w");
  io_outtrj2(outputfile,numatom,1,crd_seled);
  fclose(outputfile);

  dihed_seled=CD(crd_seled,flagKOP);
  for (i=0;i<numsdihed;++i) printf("%12.4lf\n",dihed_seled[sdihednum[i]]);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:%s \n",progname);
  printf("-s skip lines  \n");
  printf("-h help  \n");
  printf("numstep numsdihed sdihednum sdihed[-180-180] trjfilename parmtopfilename outputfilename  \n");
}



