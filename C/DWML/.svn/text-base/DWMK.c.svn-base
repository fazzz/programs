
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "RAND.h"
#include "BOXMULL.h"

#define ON 1
#define OFF 0

double cv(double x,double y, double x0, double y0, double ex, double ey, double k, double k2);
void cf(double x,double y,double ex,double ey,double k,double k2, double x0, double y0, double *fx, double *fy,double *fex,double *fey);

int main(int argc, char *argv[]) {
  int i,j,l,flag=OFF;
  int numstep;

  double dt,k=0.0/*0.0001*/,k2=0.0;
  double x,y,vx,vy,ex,vex,ey,vey;
  double x0,y0,vx0,vy0;
  double fx,fx_p,fy,fy_p,fex,fex_p,fey,fey_p;
  double ke,v,te;
  double mu,D,DE,DX,DY,DEX,DEY,kbT;

  char *trjfilename,*velfilename,*enefilename;
  FILE *trjfile,*velfile,*enefile;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  while((c=getopt(argc,argv,"B"))!=-1) {
    switch(c) {
    case 'B':
      flag=ON;
      break;
    default:
      printf("USAGE: %s dt numstep x0 y0 vx0 vy0 ex0 ey0 k k2 trjfilename velfilename enefilename \n",argv[0]);
      printf("USAGE: %s -B dt numstep x0 y0 D kbT ex ey k k2 trjfilename velfilename enefilename \n",argv[0]);
      exit(1);
    }
  }

  argc-=optind;
  argv+=optind;

  if (flag==OFF) {
    if (argc < 8) {
      printf("USAGE: %s dt numstep x0 y0 vx0 vy0 ex0 k k2 trjfilename velfilename enefilename \n",argv[0]);
      exit(1);
    }
    dt = atof(*argv);
    numstep = atoi(*++argv);
    x0  = atof(*++argv);
    y0  = atof(*++argv);
    vx0 = atof(*++argv);
    vy0 = atof(*++argv);
    ex=atof(*++argv);
    ey=atof(*++argv);
    k=atof(*++argv);
    k2=atof(*++argv);
    trjfilename  = *++argv;
    velfilename  = *++argv;
    enefilename  = *++argv;
    trjfile=efopen(trjfilename,"w");
    velfile=efopen(velfilename,"w");
    enefile=efopen(enefilename,"w");
  }
  else {
    if (argc < 10) {
      printf("USAGE: %s dt numstep x0 y0 mu kbT ex k k2 trjfilename velfilename enefilename \n",argv[0]);
      exit(1);
    }
    dt = atof(*argv);
    numstep = atoi(*++argv);
    x0  = atof(*++argv);
    y0  = atof(*++argv);
    mu  = atof(*++argv);
    kbT = atof(*++argv);
    D=kbT*mu;
    DE=/*0.1**/D;
    ex=atof(*++argv);
    ey=atof(*++argv);
    k=atof(*++argv);
    k2=atof(*++argv);
    trjfilename  = *++argv;
    velfilename  = *++argv;
    enefilename  = *++argv;
    trjfile=efopen(trjfilename,"w");
    velfile=efopen(velfilename,"w");
    enefile=efopen(enefilename,"w");
  }

  vex=0.0;
  vey=0.0;
  x=x0;
  y=y0;
  vx=vx0;
  vy=vy0;
  for (i=0;i<numstep;++i) {
    if (flag==OFF) {
      v=cv(x,y,x0,y0,ex,ey,k,k2);
      ke=0.5*vx*vx+0.5*vy*vy+0.5*vex*vex+0.5*vey*vey;
      te=ke+v;
      fprintf(trjfile,"%d %12.8lf %12.8lf %12.8lf %12.8lf\n",i+1,x,y,ex,ey);
      fprintf(velfile,"%d %12.8lf %12.8lf %12.8lf %12.8lf\n",i+1,vx,vy,vex,vey);
      fprintf(enefile,"%d %12.8lf %12.8lf %12.8lf\n",i+1,ke,v,te);
      cf(x,y,ex,ey,k,k2,x0,y0,&fx,&fy,&fex,&fey);
      fx_p=fx;									
      fy_p=fy;									
      fex_p=fex;
      fey_p=fey;
      x+=dt*vx+dt*dt/2.0*fx;
      y+=dt*vy+dt*dt/2.0*fy;
      ex+=dt*vex+dt*dt/2.0*fex;
      ey+=dt*vey+dt*dt/2.0*fey;
      cf(x,y,ex,ey,k,k2,x0,y0,&fx,&fy,&fex,&fey);
      vx+=dt/2.0*(fx+fx_p);
      vy+=dt/2.0*(fy+fy_p);
      vex+=dt/2.0*(fex+fex_p);
      vey+=dt/2.0*(fey+fey_p);
    }
    else {
      DX=2.0*D*dt*Box_Muller(i,0.0,1.0);
      DY=2.0*D*dt*Box_Muller(i,0.0,1.0);
      DEX=2.0*DE*dt*Box_Muller(i,0.0,1.0);
      DEY=2.0*DE*dt*Box_Muller(i,0.0,1.0);
      v=cv(x,y,x0,y0,ex,ey,k,k2);
      fprintf(trjfile,"%d %12.8lf %12.8lf %12.8lf %12.8lf\n",i+1,x,y,ex,ey);
      fprintf(velfile,"%d %12.8lf %12.8lf %12.8lf %12.8lf\n",i+1,vx,vy,vex,vey);
      fprintf(enefile,"%d %12.8lf %12.8lf %12.8lf\n",i+1,ke,v,te);
      cf(x,y,ex,ey,k,k2,x0,y0,&fx,&fy,&fex,&fey);
      vx=(D/kbT*fx*dt+DX)/dt;
      vy=(D/kbT*fy*dt+DY)/dt;
      vex=(DE/kbT*fex*dt+DEX)/dt;
      vey=(DE/kbT*fey*dt+DEY)/dt;
      x+=vx*dt;
      y+=vy*dt;
      ex+=vex*dt;
      ey+=vey*dt;
      ke=0.5*vx*vx+0.5*vy*vy+0.5*vex*vex+0.5*vey*vey;
      te=ke+v;
    }
  }
  fclose(trjfile);
  fclose(velfile);
  fclose(enefile);
  
  return 0;
}

double cv(double x,double y, double x0, double y0, double ex, double ey, double k, double k2) {
  double v;

  v=(x*x*x*x-2.0*x*x)+(y*y*y*y-2.0*y*y)+0.5*k*(x-y-ex)*(x-y-ex)/*0.5*k*(x-ex)*(x-ex)+0.5*k*(y-ey)*(y-ey)+0.5*k2*(ex-ey)*(ex-ey)*/+0.5*k2*ex*ex;

  return v;
}

void cf(double x,double y,double ex,double ey,double k,double k2, double x0, double y0, double *fx, double *fy,double *fex,double *fey) {
  *fx=-1.0*(4.0*x*x*x-4.0*x/*+k*(x-ex)*/+k*(x-y-ex));
  *fy=-1.0*(4.0*y*y*y-4.0*y/*+k*(y-ey)*/-k*(x-y-ex));
  *fex=/*1.0*k*(x-ex)-1.0*k2*(ex-ey)*/+k*(x-y-ex)-k2*ex;
  *fey=/*1.0*k*(y-ey)+1.0*k2*(ex-ey)*/0.0;
}
