
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "EF.h"
#include "PT.h"
#include "MB.h"

int main(int argc, char *argv[]) {
  int i,j,k,numselect;
  int numatom,numstep,atomp[8];
  double max1,max2,min1,min2,pi;
  double theta[2],old_theta[2];
  double atom[8][3];

  double *trj;
  char *inputfilename,*inputfilename2,*inputfilename3,*outputfilename;
  FILE *inputfile,*inputfile2,*inputfile3,*outputfile;

  if (argc < 5) {
    printf("USAGE: ./selecttrj.exe inputfilename(trj) inputfilename2(cond) inputfilename3(parmtop) outputfilename\n");
    exit(1);
  }
  inputfilename =  *++argv;
  inputfilename2 = *++argv;
  inputfilename3 = *++argv;
  outputfilename = *++argv;

  inputfile2=efopen(inputfilename2,"r");
  fscanf(inputfile2,"%d",&numstep);
  fscanf(inputfile2,"%d",&numatom);
  for (i=0;i<4;++i)
    fscanf(inputfile2,"%d",&atomp[i]);
  fscanf(inputfile2,"%lf",&max1);
  fscanf(inputfile2,"%lf",&min1);
  for (i=4;i<8;++i)
    fscanf(inputfile2,"%d",&atomp[i]);
  fscanf(inputfile2,"%lf",&max2);
  fscanf(inputfile2,"%lf",&min2);
  fclose(inputfile2);

  inputfile3=efopen(inputfilename3,"r");
  readParmtop(inputfile3);
  numatom=AP.NATOM;
  fclose(inputfile3);

  outputfile=efopen(outputfilename,"w");
  inputfile =efopen(inputfilename,"r");
  numselect=0;
  pi=acos(-1.0);
  for (i=0;i<numstep;++i) {
    trj=(double *)ecalloc(numatom*3,sizeof(double));
    scantraj(inputfile,trj,numatom);
    for (j=0;j<8;++j)
      for (k=0;k<3;++k)
	atom[j][k]=trj[(atomp[j]-1)*3+k];
    for (j=0;j<2;++j) {
      theta[j]=pick_dihed(atom[j*4],atom[j*4+1],atom[j*4+2],atom[j*4+3],0,old_theta[j]);
      if (theta[j] > pi)
	theta[j] = theta[j] -2.0*pi;
      theta[j] = theta[j]*180/pi;
      old_theta[j]=theta[j];
   }
    if (theta[0] < max1 && theta[0] >= min1 && theta[1] < max2 && theta[1] >= min2) {
      outtraj(outputfile,trj,numatom);
      ++numselect;
    }
    old_theta[0]=theta[0];
    old_theta[1]=theta[1];
    free(trj);
  }
  fprintf(outputfile,"%d\n",numselect);
  fclose(inputfile);
  fclose(outputfile);

}

int scantraj(FILE *inputfile, double *traj,int numatom){
  int i,j;
  double d;

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      if (fscanf(inputfile,"%lf",&d)== EOF){
	break;
      }
      traj[i*3+j]=d;
    }
  }

  return 1;

}

int outtraj(FILE *outputfile, double *selectedtraj,int numatom){
  int i,j;
  double d;

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      fprintf(outputfile,"%10.6lf ",selectedtraj[i*3+j]);
    }
    fprintf(outputfile,"\n ");
  }
  fprintf(outputfile,"\n ");

  return 1;
}


