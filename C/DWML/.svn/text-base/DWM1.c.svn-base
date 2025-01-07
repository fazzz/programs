
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int numstep;

  double dt;
  double x,vx;
  double x0,vx0;
  double fx,fx_p;
  double ke,v,te;

  char *outputfilename1,*outputfilename2;
  FILE *outputfile1,*outputfile2;

  if (argc < 7) {
    printf("USAGE: %s dt numstep x0 vx0 outputfilename1 outputfilename2 \n",argv[0]);
    exit(1);
  }
  dt = atof(*++argv);
  numstep = atoi(*++argv);
  x0 = atof(*++argv);
  vx0 = atof(*++argv);
  outputfilename1  = *++argv;
  outputfilename2  = *++argv;
  outputfile1=efopen(outputfilename1,"w");
  outputfile2=efopen(outputfilename2,"w");

  x=x0;
  vx=vx0;
  for (i=0;i<numstep;++i) {
    v=x*x*x*x-2.0*x*x;
    ke=0.5*vx*vx;
    te=ke+v;
    fprintf(outputfile1,"%d %12.8lf \n",i+1,x);
    fprintf(outputfile2,"%d %12.8lf %12.8lf %12.8lf\n",i+1,ke,v,te);
    fx=-1.0*(4.0*x*x*x-4.0*x);
    fx_p=fx;
    x+=dt*vx+dt*dt/2.0*fx;
    fx=-1.0*(4.0*x*x*x-4.0*x);
    vx+=dt/2.0*(fx+fx_p);
  }
  fclose(outputfile1);
  fclose(outputfile2);
  
  return 0;
}

