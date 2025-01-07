
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "RAND.h"
#include "BOXMULL.h"

double cv(double x,double y,double ex, double ey, double k, double k2);

int main(int argc, char *argv[]) {
  int i,j;
  int c;
  int numstep;

  double beta,delta;
  double dx=0.01,dex=0.01,k=0.0,k2=0.0;
  double x,y,ex,ey;
  double x_trial,y_trial,ex_trial,ey_trial,v_trial;
  double v;

  char *trjfilename,*enefilename,*logfilename;
  FILE *trjfile,*enefile,*logfile;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  if (argc < 8) {
    printf("USAGE: %s numstep x0 y0 ex0 ey0 k k2 trjfilename enefilename logfilename \n",argv[0]);
    exit(1);
  }
  numstep = atoi(*++argv);
  beta = atof(*++argv);
  x  = atof(*++argv);
  y  = atof(*++argv);
  ex = atof(*++argv);
  ey = atof(*++argv);
  k=atof(*++argv);
  k2=atof(*++argv);
  trjfilename  = *++argv;
  enefilename  = *++argv;
  logfilename  = *++argv;
  trjfile=efopen(trjfilename,"w");
  enefile=efopen(enefilename,"w");
  logfile=efopen(logfilename,"w");
  fprintf(logfile,"%s %s\n",trjfilename,enefilename);
  fprintf(logfile,"numstep=%d beta=%4.2lf x0=%4.2lf y0=%4.2lf ex0=%4.2lf ey0=%4.2lf \nk=%4.2lf k2=%4.2lf\n",numstep,beta,x,y,ex,ey,k,k2);

  v=cv(x,y,ex,ey,k,k2);
  for (i=0;i<numstep;++i) {
    x_trial=x+dx*Box_Muller(i,0.0,1.0);
    y_trial=y+dx*Box_Muller(i,0.0,1.0);
    ex_trial=ex+dex*Box_Muller(i,0.0,1.0);
    ey_trial=ey+dex*Box_Muller(i,0.0,1.0);
    v_trial=cv(x_trial,y_trial,ex_trial,ey_trial,k,k2);
    delta=v_trial-v;
    if((c=Metropolis(beta*delta))==1) {
      x=x_trial;
      y=y_trial;
      v=v_trial;
      ex=ex_trial;
      ey=ey_trial;
    }
    fprintf(trjfile,"%d %12.8lf %12.8lf %12.8lf %12.8lf\n",i+1,x,y,ex,ey);
    fprintf(enefile,"%d %12.8lf %d\n",i+1,v,c);
  }
  fclose(trjfile);
  fclose(enefile);
  fclose(logfile);
  
  return 0;
}

double cv(double x,double y,double ex, double ey, double k, double k2) {
  double v;

  v=(x*x*x*x-2.0*x*x)+(y*y*y*y-2.0*y*y)+0.5*k*(x-y-ex)*(x-y-ex)+0.5*k*(y-ey)*(y-ey)+0.5*k2*(ex-ey)*(ex-ey);

  return v;
}


