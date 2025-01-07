
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fftw3.h"

#include "IO.h"
#include "PT.h"
#include "EF.h"
#include "SPE.h"
#include "const.h"

void FTd(int numdihed, int nstep,double *dtrj, double *spe);
void FT(int natom, int nstep,double *trj, double *spe);

int main(int argc, char *argv[]) {
  int i,j,k,l,n;
  int flag,flag2;
  int numatom,numres,numdihed,numstep,num;
  int num_of_split,num_step_ini,num_step_fin;
  double c=2.999792e-2;
  double kb=1.98723e-3*4.18407*100.0;
  double deltat,pi,kbT,temp,sumspe;
  double *spe,*trj,*speave;

  char *inputfilename1,*inputfilename2,*inputfilename3,*outputfilename;
  FILE *inputfile1,*inputfile2,*inputfile3,*outputfile;
  
  if (argc < 5) {
    printf("USAGE: %s flag(c or v) flag2(c or a or d )  inputfilename1(vel or trj) inputfilename2(cond) inputfilename3(parmtop) outputfilename\n",argv[0]);
    printf("cond: deltat num_of_split num_step_ini num_step_fin temp (numdihed)\n");
    exit(1);
  }
  flag=(*++argv)[0];
  if (flag != 'c' && flag != 'v') {
    printf("flag error: must be c  or v ");
    exit(1);
  }
  flag2=(*++argv)[0];
  if (flag2 != 'c' && flag2 != 'a' && flag2 != 'd') {
    printf("flag error: must be c  or a or d ");
    exit(1);
  }
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  inputfilename3 = *++argv;
  outputfilename = *++argv;
  
  inputfile2=efopen(inputfilename2,"r");
  fscanf(inputfile2,"%lf",&deltat);
  fscanf(inputfile2,"%d",&num_of_split);
  fscanf(inputfile2,"%d",&num_step_ini);
  fscanf(inputfile2,"%d",&num_step_fin);  
  fscanf(inputfile2,"%lf",&temp);  
  if (flag2 == 'd')
    fscanf(inputfile2,"%d",&numdihed);        
  fclose(inputfile2);
  numstep=(num_step_fin-num_step_ini)/num_of_split;
  kbT=kb*temp;

  inputfile3=efopen(inputfilename3,"r");
  readParmtop(inputfile3);
  numres=AP.NRES;
  numatom=AP.NATOM;
  fclose(inputfile3);
  if (flag2=='d')
    num=numdihed;
  else if (flag2=='c')
    num=numres;
  else
    num=numatom;
  
  inputfile1=efopen(inputfilename1,"r");
  if (flag2=='d')
    io_dismissdata(inputfile1,num_step_ini-1,num);
  else
    io_dismissdata(inputfile1,num_step_ini-1,num*3);
  speave=(double *)ecalloc(sizeof(double),numstep*num);
  for (i=0;i<num_of_split;++i) {
    spe=(double *)ecalloc(sizeof(double),numstep*num);
    if (flag2=='d')
      trj=(double *)ecalloc(sizeof(double),numstep*num);
    else
      trj=(double *)ecalloc(sizeof(double),numstep*num*3);
    if (flag2=='d')	      
      io_scandtraj(inputfile1,numstep,num,trj);
    else if (flag2=='c')
      io_scancalphatraj_aw(inputfile1,numstep,numatom,numres,trj);
    else
      io_scantraj_aw(inputfile1,numstep,num,trj);
    if (flag2=='d')	      
      CSPd_decom(num,numstep,trj,spe);
    else
      CSP_decom(num,numstep,trj,spe);
    for (j=0;j<numstep;++j)
      for (k=0;k<num;++k)
	speave[j*num+k]=(i*speave[j*num+k]+spe[j*num+k])/(i+1);
    free(spe);
    free(trj);
  }
  fclose(inputfile1);
    
  for (n=0;n<num;++n){
    sumspe=0.0;
    for (i=0;i<numstep;++i)
      sumspe+=speave[i*num+n];
    /**********************************************************/
    /* if (flag2=='d')					      */
    /*   for (i=0;i<numstep;++i)			      */
    /* 	speave[i*num+n]=speave[i*num+n]/sumspe;		      */
    /* else						      */
    /*   for (i=0;i<numstep;++i)			      */
    /* 	speave[i*num+n]=speave[i*num+n]/sumspe*3.0;	      */
    /**********************************************************/
    for (i=0;i<numstep;++i)
      speave[i*num+n]=speave[i*num+n]/kbT;
  }

  outputfile=efopen(outputfilename,"w");
  pi=acos(-1.0);
  for (i=0;i<numstep;++i) {
    fprintf(outputfile,"%e ",(double)i/numstep/deltat/c);
    for (j=0;j<num;++j)
      if (flag=='v')
	fprintf(outputfile,"%e ",speave[i*num+j]);
      else
	fprintf(outputfile,"%e ",(2.0*pi*i/numstep/deltat)*(2.0*pi*i/numstep/deltat)*speave[i*num+j]);
    fprintf(outputfile,"\n ");
  }
  fclose(outputfile);
  free(speave);

  return 0;
}


