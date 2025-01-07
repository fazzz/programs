
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "EF.h"
#include "PT.h"
#include "MB.h"

int main(int argc, char *argv[]) {
  int i,numselect;
  int numatom,numstep,atomp[8];
  double max1,max2,min1,min2;
  double theta[2],old_theta[2];

  double *traj;
  char *inputfilename,*inputfilename2;
  FILE *inputfile,*inputfile2,*parmfile;
  FILE *outputfile;

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
  readParmTop(inputfile3);
  fclose(inputfile3);
  numatom=AP.NATOM;
  fclose(inputfile3);

  numselect=0;
  inputfile=efopen(inputfilename,"r");
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    trj=(double *)ecalloc(numatom*3,sizeof(double));
    scantraj(inputfile,trj,numatom);
    theta[0]=pick_dihed(atomp[0],atomp[1],atomp[2],atomp[3],0,old_theta[0]);
    theta[1]=pick_dihed(atomp[4],atomp[5],atomp[6],atomp[7],0,old_theta[1]);
    if (theta[0] < max1 && theta[0] >= min1 && theta[1] < max2 && theta[1] >= min2) {
      outtraj(outputfile,traj,numatom);
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
      coord_tag[i*3+j]=d;
    }
  }

  return 1;

}

int outtraj(FILE *outputfile, double *selectedtraj,int numatom){
  int i,j;
  double d;

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j)
      fprintf(outputfile,"%10.6lf",selectedtraj[i*3+j]);

  return 1;
}


