
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "FF.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0

int main(int argc, char *argv[]) {
  int i,j,k;

  double dihed0[10],dihed1[10],dihed2[10];
  double P,N[10],KE,vel[10],dt,inertia;
  double P_temp,N_temp;
  int numstep,numtype[4],numdihed,type[10][4];
  char *line,N_p;
  size_t len=0;

  char *inputfilename,*inputfilename2;
  char *outputfilename,*outputfilename2;
  FILE *inputfile,*inputfile2;
  FILE *outputfile,*outputfile2;
  FILE *log;

  if (argc < 5) {
    printf("USAGE: ./%s inertia dt numstep inputfilename1(parmtop) inputfilename2(cond) outputfilename(ene_di) outputfilename(dihed)\n",argv[0]);
    printf("cond numdihed dihed1 dihed0 numtype type1 ...... \n");
    exit(1);
  }
  inertia=atof(*++argv);
  dt=atof(*++argv);  
  numstep=atoi(*++argv);
  inputfilename  = *++argv;
  inputfilename2 = *++argv;
  outputfilename  = *++argv;
  outputfilename2  = *++argv;

  inputfile=efopen(inputfilename,"r");
  readParmtop(inputfile);
  fclose(inputfile);

  inputfile2=efopen(inputfilename2,"r");
  fscanf(inputfile2,"%d",&numdihed);
  for (i=0;i<numdihed;++i) {
    fscanf(inputfile2,"%lf %lf",&dihed1[i],&dihed0[i]);
    fscanf(inputfile2,"%d",&numtype[i]);
    for (j=0;j<numtype[i];++j) {
      fscanf(inputfile2,"%d",&type[i][j]);
    }
  }
  fclose(inputfile2);

  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");

  for (i=0;i<numstep;++i) {
    P=0.0;
    for (j=0;j<numdihed;++j)
      N[j]=0.0;
    KE=0;
    for (j=0;j<numdihed;++j) {
      for (k=0;k<numtype[j];++k) {
	P_temp=0.0;
	N_temp=0.0;
	ff_calc_spe_type_DIHE(&P_temp,&N_temp,dihed1[j],type[j][k]);
	P+=P_temp;
	N[j]+=N_temp;
      }
      dihed2[j] = 2.0*dihed1[j]-dihed0[j]+dt*dt*N[j]/inertia;
      vel[j]   = (dihed2[j]-dihed0[j])/2.0/dt;

      dihed0[j]=dihed1[j];
      dihed1[j]=dihed2[j];
      KE+=0.5*inertia*vel[j]*vel[j]/4.18407/100.0;
    }

    fprintf(outputfile, "%d  %e %e %e\n",i+1,P+KE,P,KE);
    fprintf(outputfile2,"%d ",i+1);
    for (j=0;j<numdihed;++j)
      fprintf(outputfile2,"%e ",dihed2[j]);
    fprintf(outputfile2,"\n");
  }

  fclose(outputfile);
  fclose(outputfile2);

  return 0;
}

