
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "PT.h"
#include "IO.h"
#include "rmsd.h"

#define ON 1
#define OFF 0

int main(int argc, char *argv[]) {
  int i,j;
  int flag='c',flagAmber=ON;
  int numatom,numres,numstep;
  double *crd_ref,*crd,rmsd;
  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *inputfilename1,*inputfilename2,*inputfilename3,*outputfilename;
  FILE *inputfile1,*inputfile2,*inputfile3,*outputfile;
  
  if (argc < 4 ) {
    printf("USAGE: %s [AC] numstep inputfilename1(trj) inputfilename2(top) inputfilename3(ref) \n",argv[0]);
    exit(1);
  }

  while((c=getopt(argc,argv,"AC"))!=-1) {
    switch(c) {
    case 'A':
      flagAmber=ON;
      break;
    case 'C':
      flag='a';
      break;
    default:
      printf("USAGE A-flag for Amber mode, C-flag for Ca mode");
      exit(1);
    }
  }
  argc-=optind;
  argv+=optind;
  numstep=atoi(*argv);
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  inputfilename3 = *++argv;
  
  inputfile2=efopen(inputfilename2,"r");
  readParmtop(inputfile2);
  fclose(inputfile2);
  numatom=AP.NATOM;
  numres=AP.NRES;
  
  inputfile3=efopen(inputfilename3,"r");
  if (flag=='a') {
    crd_ref=(double *)gcemalloc(sizeof(double)*numres*3);
    crd=(double *)gcemalloc(sizeof(double)*numres*3);
    io_scanconf_atomtype(inputfile3,"CA  ",2,numatom,crd_ref);
  }
  else {
    crd_ref=(double *)gcemalloc(sizeof(double)*numatom*3);
    crd=(double *)gcemalloc(sizeof(double)*numatom*3);
    io_scanconf(inputfile3,numatom,crd_ref,'x');
  }
  fclose(inputfile3);

  inputfile1=efopen(inputfilename1,"r");
  if (flagAmber==ON)
    getline(&line,&len,inputfile1);
  for (i=0;i<numstep;++i) {
    if (flag=='a') {
      io_scanconf_atomtype(inputfile1,"CA  ",2,numatom,crd);
    }
    else {
      io_scanconf(inputfile1,numatom,crd,'x');
    }

    rmsd=rmsd_qcp(crd,crd_ref,numres);
    printf("%d %e\n",i+1,rmsd);
  }
  fclose(inputfile1);
  
  return 0;
}

