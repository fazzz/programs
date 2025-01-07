
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PT.h"
#include "IO.h"
#include "BAR.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j;

  double *traj1,*traj2;

  int n1,n2,natom;
  double C,**ebU1,**ebU2,dF;
  
  char *inputfilename1,*inputfilename2,*inputfilename3,*outputfilename;
  FILE *inputfile1,*inputfile2,*inputfile3,*outputfile;
  
  if (argc < 4) {
    printf("USAGE: %s inputfilename1(trj) inputfilename2(parm) inputfilename3(cond) outputfilename\n",argv[0]);
    exit(1);
  }
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  inputfilename3 = *++argv;
  outputfilename = *++argv;
  
  inputfile3=efopen(inputfilename3,"r");
  fscanf(inputfile3,"%d",&n1);
  fscanf(inputfile3,"%d",&n2);
  fclose(inputfile3);

  inputfile4=efopen(inputfilename4,"r");
  readParmtop(inputfile4);
  fclose(inputfile4);
  natom=AP.NATOM;

  traj1=(double *)gcemalloc(sizeof(double)*n1*natom*3);
  traj2=(double *)gcemalloc(sizeof(double)*n2*natom*3);

  ebU1[0]=(double *)gcemalloc(sizeof(double)*n1);
  ebU1[1]=(double *)gcemalloc(sizeof(double)*n2);
  ebU2[0]=(double *)gcemalloc(sizeof(double)*n1);
  ebU2[1]=(double *)gcemalloc(sizeof(double)*n2);
  
  inputfile1=efopen(inputfilename1,"r");
  io_scantraj(inputfile1,n1,natom,traj1);
  fclose(inputfile1);

  inputfile2=efopen(inputfilename2,"r");
  io_scantraj(inputfile2,n2,natom,traj2);
  fclose(inputfile2);

  BAR_CCF(n1,n2,traj1,traj2,ebU1,ebU2);
  BAR_ite(n1,n2,C,ebU1,ebU2,dF);
  
  outputfile=efopen(outputfilename,"w");
  fprintf(outputfile,"C=%lf\n",C);
  fprintf(outputfile,"dF=%lf\n",dF);
  fclose(outputfile);
  
  return 0;
}

