
#include <stdio.h>

#include <stdlib.h>

#include <math.h>

#include "FourierTransform.h"
#include "ParmTop.h"

int scandtraj(FILE *inputfile, double *dtraj, int numdihed, int numstep);

int main(int argc, char *argv[]) {
  int i,j,k;
  int numdihed,numstep;
  double deltat,pi;
  double *dtraj,*Ftraj_R,*Ftraj_I,*PowerSpectral;
  char *inputfilename,*inputfilename2,*outputfilename;
  FILE *inputfile,*inputfile2, *outputfile;

  if (argc < 3) {
    printf("USAGE: ./calc_dihed_spectral.exe  inputfilename(dihed) inputfilename(bond length, bond angle......)  outputfilename(coordinate)\n");
    exit(1);
  }
  inputfilename =  *++argv;
  inputfilename2 = *++argv;
  outputfilename = *++argv;

  if((inputfile2=fopen(inputfilename2,"r"))==NULL)  {
    printf("There is not %s\n",inputfilename2);
    exit(1);
  }
  fscanf(inputfile2,"%d",&numstep);
  fscanf(inputfile2,"%d",&numdihed);
  fscanf(inputfile2,"%lf",&deltat);
  fclose(inputfile2);
  
  if((inputfile=fopen(inputfilename,"r"))==NULL) {
    printf("There is not %s\n",inputfilename);
    exit(1);
  }
  dtraj = malloc(sizeof(double)*numstep*numdihed);
  scandtraj(inputfile,dtraj,numdihed,numstep);
  fclose(inputfile);

  
  if((outputfile=fopen(outputfilename,"w"))==NULL) {
    printf("There is not %s\n",outputfilename);
    exit(1);
  }
  
  fclose(outputfile);
  free(dtraj);
  free(Ftraj_R);
  free(Ftraj_I);
  free(PowerSpectral);
  
  return 0;
}
 
int scandtraj(FILE *inputfile, double *dtraj, int numdihed, int numstep) {
   int i,j;
   
   for (i=0;i<numstep;++i)  {
     for (j=0;j<numdihed;++j)  {
       fscanf(inputfile,"%lf",&dtraj[i*numdihed+j]);
     }
   }

   return 0;
 }
 
