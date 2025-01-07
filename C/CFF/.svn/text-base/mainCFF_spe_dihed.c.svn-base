
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
  int i;

  double P,N,KE,angle,vel,dt,inertia,angle_p,angle_pp;
  int numstep,type;
  char *line,N_p;
  size_t len=0;

  char *inputfilename;
  char *outputfilename,*outputfilename2;
  FILE *inputfile;
  FILE *outputfile,*outputfile2;
  FILE *log;

  if (argc < 5) {
    printf("USAGE: ./%s angle0 angle1 inertia type dt numstep inputfilename1(parmtop) outputfilename(ene_di)\n",argv[0]);
    exit(1);
  }
  angle_p=atof(*++argv);
  angle_pp=atof(*++argv);
  inertia=atof(*++argv);
  type=atoi(*++argv);
  dt=atof(*++argv);  
  numstep=atoi(*++argv);
  inputfilename  = *++argv;
  outputfilename  = *++argv;
  outputfilename2  = *++argv;

  inputfile=efopen(inputfilename,"r");
  readParmtop(inputfile);
  fclose(inputfile);

  if (type >= AP.NUMANG)
    printf("error! nothing such type......\n");

  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");

  P=0.0;
  N=0.0;
  for (i=0;i<numstep;++i) {
    ff_calc_spe_type_DIHE(&P,&N,angle_p,type);

    angle = 2.0*angle_p-angle_pp+dt*dt*N/inertia;
    vel   = (angle-angle_pp)/2.0/dt;

    angle_pp=angle_p;
    angle_p=angle;
    KE=0.5*inertia*vel*vel/4.18407/100.0;

    fprintf(outputfile, "%d  %e %e %e\n",i+1,P+KE,P,KE);
    fprintf(outputfile2,"%d  %e \n",i+1,angle);
  }

  fclose(outputfile);
  fclose(outputfile2);

  log=efopen("log_spe_dihed.txt","w");
  fprintf(log,"type=%d V=%lf n=%lf theta=%lf\n",type,AP.PK[type],AP.PN[type],AP.PHASE[type]);
  fclose(log);


  return 0;
}

