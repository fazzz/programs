

#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

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
  
  char *inputfilename,*outputfilename,*parmtopfilename;
  FILE *inputfile,*outputfile,*parmtopfile;
  
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
  
  if (argc < ) {
    USAGE(progname);
    exit(1);
  }
  inputfilename=*++argv;
  parmtopfilename=*++argv;
  outputfilename=*++argv;
  
  inputfile=efopen(inputfilename,"r");
  for (i=0;i<sl;++i) getline(&line,&len,inifile);
  io_scantrj2(inifile,numatom,numstep,path);
  fclose(inputfile);

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);

  for (i=0;i<numstep;++i) {
    for (j=0;j<numatom*3;++j)  {
      crd[j]=path[i*numatom*3+j];
      CD(crd,dihed,flag);
    }
  }

  sort_by_dihedral();
  
  outputfile=efopen(outputfilename,"w");
  io_outtrj2(outputfile,numatom,numpoint,path);  
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:%s \n",progname);
  printf("-h help  \n");
  printf("inputfilename outputfilename  \n");
}



