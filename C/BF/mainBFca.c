
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PT.h"
#include "BF.h"
#include "IO.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j,k,dummy;
  int numatom,numres,numstep,rstflag;
  double *mass,*crd_ref,*crd;
  double rmsd=0.0;
  char *inputfilename,*inputfilename2,*inputfilename3,*outputfilename;
  FILE *inputfile,*inputfile2,*inputfile3,*outputfile;
  char *line;
  size_t len=0;

  if (argc < 5) {
    printf("USAGE:./carmsd.exe numstep rstflag input(traj) input(ref crdfile ) input(ParmTop) outputname(rmsd)\n");
    exit(1);
  }
  numstep=atoi(*++argv);
  rstflag=atoi(*++argv);
  inputfilename   =  *++argv;
  inputfilename2  =  *++argv;
  inputfilename3  =  *++argv;
  outputfilename  =  *++argv;

  inputfile3=efopen(inputfilename3,"r");
  readParmtop(inputfile3);
  fclose(inputfile3);
  numatom = AP.NATOM;
  numres=AP.NRES;
  mass=(double *)gcemalloc(sizeof(double)*numres);
  crd=(double *)gcemalloc(sizeof(double)*numres*3);
  crd_ref=(double *)gcemalloc(sizeof(double)*numres*3);
  j=0;
  for (i=0;i<numatom;++i) {
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
     mass[j] = AP.AMASS[i];     ++j;
    }
  }

  inputfile2=efopen(inputfilename2,"r");
  fscanf(inputfile2,"%d",&dummy);
  if (rstflag==1)
    fscanf(inputfile2,"%d",&dummy);
  io_scanconf_atomtype(inputfile2,"CA",2,numatom,crd_ref);
  fclose(inputfile2);

  outputfile=efopen(outputfilename,"w");
  inputfile=efopen(inputfilename,"r");

  getline(&line,&len,inputfile);
  for (i=0;i<numstep;++i) {
    io_scanconf_atomtype(inputfile,"CA",2,numatom,crd);
    rmsd = bestfit_crd(crd_ref,crd,mass,numres);
    fprintf(outputfile,"%d %e\n",i,rmsd);
  }
  fclose(inputfile);
  fclose(outputfile);

  return 0;
}


