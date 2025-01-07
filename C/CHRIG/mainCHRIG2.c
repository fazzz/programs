#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PT.h"
#include "FF.h"
#include "EF.h"
#include "EF.h"
#include "QUA.h"
#include "TOPO.h"
#include "MB.h"

#include "netcdf_mine.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,d;
  int numstep,numatom;
  double dt=0.001;
  double U=418.4070;
  double q[4],axis[3],length;
  double delta_dihed;
  double total,*hist,*pmf,*p_now,*dihed,pmf_min;
  double atom1[3],atom2[3],atom3[3],atom4[3];

  double *crd,crd_now[7][4],crd_rot[7][4];
  double *ene;
  struct potential e;
  struct force f;

  double pi;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*outputfilename,*outputfilename2,*parmfilename;
  FILE *inputfile,*outputfile,*outputfile2,*parmfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"h",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  numstep = atoi(*argv);
  //  inputfilename = *++argv;
  //  parmfilename = *++argv;
  outputfilename = *++argv;
  outputfilename2 = *++argv;

  pi=acos(-1.0);

  inputfilename="/home/yamamori/calspa/input/BT.crd";
  parmfilename="/home/yamamori/calspa/input/BT.top";

  delta_dihed=2.0*acos(-1.0)/numstep;

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
 
  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  inputfile=efopen(inputfilename,"r");
  getline(&line,&len,inputfile);
  fscanf(inputfile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfile,"%lf",&crd[i*3+j]);
  fclose(inputfile);

  for(i=0;i<3;++i) axis[i] = crd[2*3+i]-crd[1*3+i];
  length = 0.0;
  for(i=0;i<3;++i) length += axis[i]*axis[i]; 
  length = sqrt(length);
  for(i=0;i<3;++i) axis[i] = axis[i]/length;

  q[0]=cos(delta_dihed*0.5);
  for(i=0;i<3;++i) q[i+1] = axis[i]*sin(delta_dihed*0.5);

  ff_set_calcffandforce(&e,&f);
  for (i=0;i<7;++i) crd_now[i][0]=0.0;

  p_now=(double *)gcemalloc(sizeof(double)*numstep);
  dihed=(double *)gcemalloc(sizeof(double)*numstep);

  outputfile=efopen(outputfilename,"w");
  total=0;
  for (i=0;i<numstep;++i) {
    ff_calcffandforce(crd,numatom,&e,&f);
    for (j=0;j<3;++j) {
      atom1[j]=crd[j];
      atom2[j]=crd[3+j];
      atom3[j]=crd[6+j];
      atom4[j]=crd[9+j];
    }
    //    dihed=dih(atom1,atom2,atom3,atom4);
    dihed[i]=pick_dihed(atom1,atom2,atom3,atom4,0,0);
    if (dihed[i]<-1.0*pi) dihed[i]+=2.0*pi;
    if (dihed[i]>pi) dihed[i]-=2.0*pi;

    for (k=0;k<3;++k) {
      crd_now[0][k+1]  = crd[2*3+k]-crd[1*3+k];
      crd_now[1][k+1]  = crd[3*3+k]-crd[1*3+k];
      crd_now[2][k+1]  = crd[9*3+k]-crd[1*3+k];
      crd_now[3][k+1] = crd[10*3+k]-crd[1*3+k];
      crd_now[4][k+1] = crd[11*3+k]-crd[1*3+k];
      crd_now[5][k+1] = crd[12*3+k]-crd[1*3+k];
      crd_now[6][k+1] = crd[13*3+k]-crd[1*3+k];
    }
    for (j=0;j<7;++j)
      qua_rot(crd_now[j],q,crd_rot[j]);
    for (k=0;k<3;++k) {
      crd[2*3+k]  = crd_rot[0][k+1]+crd[1*3+k];
      crd[3*3+k]  = crd_rot[1][k+1]+crd[1*3+k];
      crd[9*3+k]  = crd_rot[2][k+1]+crd[1*3+k];
      crd[10*3+k] = crd_rot[3][k+1]+crd[1*3+k];
      crd[11*3+k] = crd_rot[4][k+1]+crd[1*3+k];
      crd[12*3+k] = crd_rot[5][k+1]+crd[1*3+k];
      crd[13*3+k] = crd_rot[6][k+1]+crd[1*3+k];
    }

    fprintf(outputfile,"%8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf\n",dihed[i]*180/pi,e.p_t,e.p_e_t,e.p_e_14_t,e.p_LJ_t,e.p_LJ_14_t,e.p_d_t,e.p_a_t,e.p_b_t);
    total+=exp(-1.0*e.p_t);
    p_now[i]=e.p_t;
  }
  fclose(outputfile);

  hist=(double *)gcemalloc(sizeof(double)*numstep);
  for (i=0;i<numstep;++i) hist[i]=exp(-1.0*p_now[i])/total;
  
  pmf=(double *)gcemalloc(sizeof(double)*numstep);
  pmf_min=0.0;
  for (i=0;i<numstep;++i) {
    if ( pmf_min == 0.0 && hist[i] > 0.0 ) pmf_min=hist[i];
    if ( pmf_min > pmf[i] && hist[i] > 0.0 ) pmf_min=hist[i];
  }

  for (i=0;i<numstep;++i) if (hist[i]!=0.0) pmf[i]=-log(hist[i])+log(pmf_min);

  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<numstep;++i) fprintf(outputfile2,"%8.3lf %8.3lf \n",dihed[i]*180/pi,pmf[i]);
  fclose(outputfile2);
    
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] outputfilename \n",progname);
}
