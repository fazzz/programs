// Principal Component Analysis (Dihedral Angle Version)
// 03_10

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>

#include "const.h"
#include "dPCA.h"
#include "ParmTop.h"

#define MAXFILE 10;

double *traj;
//double *dihed_traj;
double *dihed_ave;


void scantraj(FILE *inputfile, /*double *traj,*/int time, int numatom);

int main(int argc, char *argv[])
{
  int i,j,k;
  int minnumoption=4;
  char *inputfilename0,*inputfilename1,*inputfilename2, *outputfilename;
  char *USAGE="USAGE: ./calc_dPCA condfile trajfile parmfile outputfile\n";
  char *ERROR="does not exist in this directry\n";
  FILE *inputfile0,*inputfile1,*inputfile2,*outputfile;

  int numatom,numdihed,time,num;
  int atom_dihed_pair[MAXNUMDIHED][4];
  double covR[MAXNUMDIHED*MAXNUMDIHED];
  double w[MAXNUMDIHED];
  double sum,cont;


  if (argc < minnumoption){
    printf("%s",USAGE);
    exit(1);
  }

  inputfilename0 = *++argv;
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  outputfilename = *++argv;

  if ((inputfile0=fopen(inputfilename0,"r"))==NULL){
    printf("%s %s",inputfilename0, ERROR);
    exit(1);
  }
  fscanf(inputfile0,"%d",&time);
  fscanf(inputfile0,"%d",&numatom);
  fscanf(inputfile0,"%d",&numdihed);
  for (i=0;i<numdihed;++i){
    for (j=0;j<4;++j){
      fscanf(inputfile0,"%d",&num);
      atom_dihed_pair[i][j]=num-1;
    }
  }
  fclose(inputfile0);

  if ((inputfile1=fopen(inputfilename1,"r"))==NULL){
    printf("%s %s",inputfilename1, ERROR);
    exit(1);
  }
  scantraj(inputfile1,/*traj,*/time,numatom);
  fclose(inputfile1);

  if ((inputfile2=fopen(inputfilename2,"r"))==NULL){
    printf("%s %s",inputfilename2, ERROR);
    exit(1);
  }
  readParmtop(inputfile2);

  //  cov=malloc(sizeof(double)*numatom*3*numatom*3);
  dihed_ave=malloc(sizeof(double)*numdihed);
  pick_dihed(/*dihed_traj,*/traj,numatom,numdihed,atom_dihed_pair,time);
  calc_dPCA(/*dihed_traj,*//*cov,*/time,numdihed,covR,w);
  proj_dPCA(/*dihed_traj,*//*traj_trans,*/dihed_ave,time,numdihed,covR,w);

  if ((outputfile=fopen("eigen_value","w"))==NULL){
    printf("eigen_value %s",ERROR);
    exit(1);
  }

  sum=0.0;
  for (i=0;i<numdihed;++i)
    sum+=w[i];

  cont=0.0;
  for (i=numdihed-1;i>=0;--i){
    cont+=w[i];
    fprintf(outputfile,"%4d -th  %12.8lf %12.8lf \n",numatom*3-i, w[i], cont/sum);
  }

  for (j=0;j<numdihed;++j){
    fprintf(outputfile,"%8.3lf",covR[numdihed-1+j*numdihed]);
  }
 fprintf(outputfile,"\n");
 for (j=0;j<numdihed;++j){
   fprintf(outputfile,"%8.3lf",covR[numdihed-2+j*numdihed]);
  }
 fprintf(outputfile,"\n");
 fclose(outputfile);

  if ((outputfile=fopen(outputfilename,"w"))==NULL){
    printf("%s %s",outputfilename,ERROR);
    exit(1);
  }

  if ((outputfile=fopen(outputfilename,"w"))==NULL){
    printf("%s %s",outputfilename,ERROR);
    exit(1);
  }

  for (i=0;i<time;++i) {
    fprintf(outputfile,"%lf %lf\n",dihed_traj_trans[i*2/*numatom*3*/],dihed_traj_trans[i*2/*numatom*3*/+1]);
  }

  fclose(outputfile);

  return 0;

}

void scantraj(FILE *inputfile, /*double *traj,*/int time, int numatom)
{
  int i,j,k;

  traj=malloc(sizeof(double)*time*numatom*3);

  for (i=0;i<time;++i)
    for (j=0;j<numatom;++j)
      for (k=0;k<3;++k)
	fscanf(inputfile,"%lf", &traj[i*numatom*3+j*3+k]);

}

