
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fftw3.h"

#include "IO.h"
#include "PT.h"
#include "EF.h"
//#include "SPE.h"
//#include "const.h"

#define kbkcl 1.98723e-3
#define kbuap 1.98723e-3*4.18407*100.0

void FTd(int numdihed, int nstep,double *dtrj, double *spe);
void FT(int natom, int nstep,double *trj, double *spe);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int flag,flag2;
  int numatom,numres,numdihed,numstep;
  int num_of_split,num_step_ini,num_step_fin;
  double c=2.999792e-2;
  double kb=1.98723e-3*4.18407*100.0;
  double deltat,pi,kbT,temp,sumspe,dof;
  double *spe,*trj,*speave;

  char *inputfilename1,*inputfilename2,*inputfilename3,*outputfilename;
  FILE *inputfile1,*inputfile2,*inputfile3,*outputfile;
  FILE *log;
  
  if (argc < 5) {
    printf("USAGE: ./CSPE.exe flag(c or v) flag2(c or a or d )  inputfilename1(vel or trj) inputfilename2(cond) inputfilename3(parmtop) outputfilename\n");
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
  if (flag=='v')
    kbT=kbkcl*temp;
  else
    kbT=kb*temp;

  inputfile3=efopen(inputfilename3,"r");
  readParmtop(inputfile3);
  numres=AP.NRES;
  numatom=AP.NATOM;
  fclose(inputfile3);
  
  inputfile1=efopen(inputfilename1,"r");
  if (flag2=='d')
    io_dismissdata(inputfile1,num_step_ini-1,numdihed);
  else if (flag2=='c')
    io_dismisstrj(inputfile1,num_step_ini-1,/*numres*/numatom);
  else if (flag2=='a')
    io_dismisstrj(inputfile1,num_step_ini-1,numatom);
  speave=(double *)ecalloc(sizeof(double),numstep);
  for (i=0;i<num_of_split;++i) {
    spe=(double *)ecalloc(sizeof(double),numstep);
    if (flag2=='d') {
      trj=(double *)ecalloc(sizeof(double),numstep*numdihed);
      io_scandtraj(inputfile1,numstep,numdihed,trj);
      FTd(numdihed,numstep,trj,spe);
    }
    else if (flag2=='a') {
      trj=(double *)ecalloc(sizeof(double),numstep*numatom*3);
      io_scantraj_aw(inputfile1,numstep,numatom,trj);
      FT(numatom,numstep,trj,spe);
    }
    else {
      trj=(double *)ecalloc(sizeof(double),numstep*numres*3);
      io_scancalphatraj_aw(inputfile1,numstep,numatom,numres,trj);
      FT(numres,numstep,trj,spe);
    }
    for (j=0;j<numstep;++j)
      speave[j]+=spe[j];
    free(spe);
    free(trj);
  }
  for (i=0;i<numstep;++i)
    speave[i]=speave[i]/num_of_split;
  fclose(inputfile1);

  pi=acos(-1.0);
    
  /***************************/
  /* sumspe=0.0;	     */
  /* for (i=0;i<numstep;++i) */
  /*   sumspe+=speave[i];    */
  /***************************/
  /* dof=0.0;								     */
  /* if (flag!='v')							     */
  /*   for (i=0;i<numstep/2;++i)					     */
  /*     dof+=(2.0*pi*i/numstep/deltat)*(2.0*pi*i/numstep/deltat)*speave[i]; */
  /* 									     */
  /* dof=sumspe/kbkcl/temp;						     */
  /***************************************************************************/

  dof=0.0;
  log=efopen("log_CSPE.txt","w");
  for (i=0;i<(int)(numstep/2);++i) {
    if (flag=='v') {
      dof+=speave[i]/kbT*2;
    }
    else {
      dof+=(2.0*pi*i/numstep/deltat)*(2.0*pi*i/numstep/deltat)*speave[i]/kbT*2;
    }
  } 
  fprintf(log,"dof=%lf\n",dof);
  fclose(log);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    if (flag=='v') {
      fprintf(outputfile,"%e %e \n",(double)i/numstep/deltat/c,speave[i]/kbT*2);
    }
    else {
      fprintf(outputfile,"%e %e \n",(double)i/numstep/deltat/c,(2.0*pi*i/numstep/deltat)*(2.0*pi*i/numstep/deltat)*speave[i]/kbT*2);
    }  
  } 
  fclose(outputfile);
  free(speave);

  
  return 0;
}

void FT(int natom, int nstep,double *trj, double *spe){
  int i,j,k,l;
  fftw_complex *in,*out;
  fftw_plan p;
  FILE *outputfile;

  double sumspe,K=0.0,sum=0.0;

  for (i=0;i<nstep;++i)
    spe[i]=0.0;

  for (i=0;i<natom;++i) {
    for (j=0;j<3;++j) {
      in  = fftw_malloc(sizeof(fftw_complex)*nstep);
      out = fftw_malloc(sizeof(fftw_complex)*nstep);
      p   = fftw_plan_dft_1d(nstep,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
      
      for (k=0;k<nstep;++k) {
	in[k][0]=trj[k*natom*3+i*3+j];
	in[k][1]=0.0;
	out[k][0]=0.0;
	out[k][1]=0.0;
      }
      fftw_execute(p);
      
      for (k=0;k<nstep;++k)
	spe[k] += out[k][0]/nstep*out[k][0]/nstep+out[k][1]/nstep*out[k][1]/nstep;
            
      fftw_destroy_plan(p);
      fftw_free(in); fftw_free(out);
    }
  }    

  for (k=0;k<nstep*natom*3;++k)
    K+=trj[k]*trj[k]/nstep;
  for (k=0;k<nstep;++k)
    sum += spe[k];


}

void FTd(int numdihed, int nstep,double *dtrj, double *spe){
  int i,j,k;
  fftw_complex *in,*out;
  fftw_plan p;
  FILE *outputfile;

  for (i=0;i<nstep;++i)
    spe[i]=0.0;

  for (i=0;i<numdihed;++i) {
    in  = fftw_malloc(sizeof(fftw_complex)*nstep);
    out = fftw_malloc(sizeof(fftw_complex)*nstep);
    p   = fftw_plan_dft_1d(nstep,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

    for (j=0;j<nstep;++j)
      for (k=0;k<2;++k)
	in[j][k]=0.0;

    for (j=0;j<nstep;++j)
	in[j][0]=dtrj[j*numdihed+i];
    
    fftw_execute(p);
    
    for (j=0;j<nstep;++j)
      spe[j] += out[j][0]*out[j][0]+out[j][1]*out[j][1];
    
    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
  }  
}
