
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SHELL.h"
#include "EF.h"


int mkdir(char *name) {
  FILE *shell;
  char *filename="/home/yamamori/mybin/mkdirsh";

  shell=efopen(filename,"w");
  fprintf(shell,"#!/bin/sh\n if [ ! -f %s ]; then \n     mkdir %s \n fi\n",name,name);
  fclose(shell);
  system("chmod +x /home/yamamori/mybin/mkdirsh");
  system("/home/yamamori/mybin/mkdirsh");
  return 1;
}




int mkmdin(int ensflag, int rstflag, int ftflag,int ftflag2, int mad, int pepcaflag,int odf, int ndf, int vof, double tp, double qi, double dt, int tl,  int os, int ose, char *name) {
  FILE *mdin;

  mdin=efopen(name,"w");
  fprintf(mdin," %d %d %lf 0.0 %lf %d 1 1 %d %d 0 0 0.0 %lf 1.0 %d %d %d 0 0 0.0 %d %d %d %d",ensflag, rstflag,tp,dt,tl,os,ose,qi,odf,ndf,vof,ftflag,ftflag2,mad,pepcaflag);
  fclose(mdin);
  return 1;
}

int makeqsub(char *inpmd,char *crd, char *vel,char *top, char *clt, char *rst, char *out, char *pn, char *dir, char *pepca, char *qsubfilename){
  FILE *qsubsh;
  char qsubfilenamefullpath[100];

  sprintf(qsubfilenamefullpath,"%s/%s",dir,qsubfilename);
  qsubsh=efopen(qsubfilenamefullpath,"w");

  fprintf(qsubsh,"#!/bin/sh\n");
  fprintf(qsubsh,"#PBS -N %s_remd\n",pn);
  fprintf(qsubsh,"#PBS -q serial\n");
  fprintf(qsubsh,"#PBS -l nodes=1\n");
  fprintf(qsubsh,"#PBS -l ncpus=1\n");
  fprintf(qsubsh,"#PBS -e %s/%s_remd.e\n",dir,pn);
  fprintf(qsubsh,"#PBS -o %s/%s_remd.o\n",dir,pn);
  fprintf(qsubsh,"cd %s\n",dir);
  fprintf(qsubsh,"/home/yamamori/mybin/MD.exe -m %s -c %s -d %s -p %s -r %s ",inpmd,crd,vel,top,clt);  
  fprintf(qsubsh,"-o %s.out -x %s.trj -v %s.vel -j %s.rst -l %s.rve -k %s.rsa -e %s\n",out,out,out,rst,rst,out,pepca);

  fclose(qsubsh);
  return 1;
}

int wait(char *name, int time) {
  FILE *SHELL;
  char *filename="/home/yamamori/mybin/wait.sh";

  SHELL=efopen(filename,"w");
  fprintf(SHELL,"#!/bin/sh\n sleep %d\nflag=`qstat | grep -G %s |wc -l`\nwhile [ ${flag} != 0 ]; do\nsleep %d\nflag=`qstat | grep -G %s |wc -l`\ndone",time,name,time,name);
  fclose(SHELL);
  
  system("chmod +x /home/yamamori/mybin/wait.sh");
  system(filename);

  return 1;
}

int makeqsub_serial_Amber(char *inpmd,char *crd,char *ref, char *vel,char *top, char *rst, char *trj,char *out, char *pn, char *dir,  char *qsubfilename, char *simtype){
  FILE *qsubsh;
  char qsubfilenamefullpath[100];

  sprintf(qsubfilenamefullpath,"%s/%s",dir,qsubfilename);
  qsubsh=efopen(qsubfilenamefullpath,"w");

  fprintf(qsubsh,"#!/bin/sh\n");
  fprintf(qsubsh,"#PBS -N %s\n",simtype);
  fprintf(qsubsh,"#PBS -q serial\n");
  fprintf(qsubsh,"#PBS -l nodes=1\n");
  fprintf(qsubsh,"#PBS -l ncpus=1\n");
  fprintf(qsubsh,"#PBS -e %s/%s.e\n",dir,simtype);
  fprintf(qsubsh,"#PBS -o %s/%s.o\n",dir,simtype);
  fprintf(qsubsh,"cd %s\n",dir);
  fprintf(qsubsh,"export AMBERHOME=/home/software/amber/amber9\n"); 
  fprintf(qsubsh,"/home/software/amber/amber9/exe/sander -O -i %s -o %s -c %s -p %s -r %s ",inpmd,out,crd,top,rst);
  if(strncmp(simtype,"min",3)==0)
    fprintf(qsubsh," -ref %s\n",ref);
  if(strncmp(simtype,"rex",3)==0)
    fprintf(qsubsh," -ref %s -x %s  -v %s\n",ref,trj,vel);
  if(strncmp(simtype,"sam",3)==0)
    fprintf(qsubsh," -x %s  -v %s\n",trj,vel);

  fclose(qsubsh);
  return 1;
}

void makeqsubbase(char *dircommand,char *dirlog,char *dir, 
		  char *name, int numnode, 
		  char *log, char *qsubfilename){
  int i;
  char *qsubsh;
  char qsubfilenamefullpath[100];

  sprintf(qsubfilenamefullpath,"%s/qsub_%s.sh",dircommand,name);
  qsubsh=efopen(qsubfilenamefullpath,"w");

  fprintf(qsubsh,"#!/bin/sh\n");
  fprintf(qsubsh,"#$ -S /bin/sh\n");
  fprintf(qsubsh,"#$ -cwd\n");
  fprintf(qsubsh,"#$ -N %s\n",name);
  fprintf(qsubsh,"#$ -q all.q@ajisai01\n");
  fprintf(qsubsh,"#$ -q all.q@ajisai0%d\n",numnode);
  fprintf(qsubsh,"#$ -e %s/%s.e\n",dirlog,name);
  fprintf(qsubsh,"#$ -o %s/%s.o\n",dirlog,name);
  fprintf(qsubsh,"cd %s\n",dir);
  fprintf(qsubsh,"LD_LIBRARY_PATH=/opt/intel/Compiler/11.0/083/mkl/lib/em64t:/opt/intel/Compiler/11.0/083/lib/intel64:/opt/intel/Compiler/11.0/083/lib/intel64:/opt/intel/impi/3.2.1.009/lib64:/opt/intel/Compiler/11.0/083/mkl/lib/em64t:/opt/intel/Compiler/11.0/083/lib/intel64:/opt/intel/Compiler/11.0/083/lib/intel6\n");
  fclose(qsubsh);
}
