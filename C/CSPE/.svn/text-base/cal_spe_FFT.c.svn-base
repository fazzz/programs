/* spectral calculation program */
/*            05_2010           */
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fftw3.h"

#include "ParmTop.h"

#define MAXSPLIT 20

double c=2.999792e-2;
double *abF;

void FourierTransformTraj(int natom, int nstep);
void FourierTransformdTraj(int numdihed, int nstep);
int scantraj(FILE *inputfile,int numatom, int N);
int scan_calpha_traj(FILE *inputfile,int numatom, int numres, int N);
int scandtraj(FILE *inputfile,int numdihed, int N);

double *traj;

int main(int argc, char *argv[]) {
  int i,j,k,DOF,ns;
  int N,M,num_of_split,numdihed,num_of_split_ini,num_of_split_fin;
  int flag=0,flag2=0;
  double deltat=1.0,pi,ab=0.0,ABS,temp;
  double *abF_ave;
  double *PowerSpectral;

  char *inputfilename,*inputfilename2,*inputfilename3,*outputfilename;
  FILE *inputfile,*inputfile2, *outputfile, *parmfile;
  char *line;
  size_t len=0;

  if (argc < 5) {
    printf("USAGE: ./calcspectral.exe flag(coordinate or velocity) flag(calpha or all atom or dihed)  inputfilename(vel or trj) inputfilename2(cond)  inputfilename3(parmtop) outputfilename\n");
    exit(1);
  }
  flag=(*++argv)[0];
  if (flag != 'c' && flag != 'v') {
    printf("flag error: must be c or v");
    exit(1);
  }
  flag2=(*++argv)[0];
  if (flag2 != 'c' && flag2 != 'a'&& flag2 != 'd') {
    printf("flag2 error: must be c or a or d");
    exit(1);
  }
  inputfilename =  *++argv;
  inputfilename2 = *++argv;
  inputfilename3 = *++argv;
  outputfilename = *++argv;

  if((inputfile2=fopen(inputfilename2,"r"))==NULL)  {
    printf("There is not %s\n",inputfilename2);
    exit(1);
  }
  fscanf(inputfile2,"%d",&M);
  fscanf(inputfile2,"%d",&N);
  fscanf(inputfile2,"%lf",&deltat);
  fscanf(inputfile2,"%lf",&temp);
  fscanf(inputfile2,"%d",&DOF);
  fscanf(inputfile2,"%d",&num_of_split);
  fscanf(inputfile2,"%d",&num_of_split_ini);
  fscanf(inputfile2,"%d",&num_of_split_fin);
  fclose(inputfile2);
  numdihed=DOF;

  if((PowerSpectral=(double *)malloc(sizeof(double)*num_of_split*M))==NULL) {
    printf("allocation error!!\n");
    exit(1);
  }

  if((parmfile=fopen(inputfilename3,"r"))==NULL) {
    printf("There is not %s\n",inputfilename3);
    exit(1);
  }  
  readParmtop(parmfile);
  fclose(parmfile);
  if (flag2=='c') {
    N=AP.NRES;
  }
  else {
    N=AP.NATOM;
  }

  if((inputfile=fopen(inputfilename,"r"))==NULL) {
    printf("There is not %s\n",inputfilename);
    exit(1);
  }
  getline(&line,&len,inputfile);
  if((outputfile=fopen(outputfilename,"w"))==NULL) {
    printf("There is not %s\n",outputfilename);
    exit(1);
  }


  abF_ave = (double *)malloc(sizeof(double)*M);
  for (i=0;i<M;++i) {
    abF_ave[i]=0.0;
  }

  for (ns=0;ns<num_of_split;++ns) {
    if (flag2=='a')
      ab=scantraj(inputfile,N,M);
    else if (flag2=='c')
      ab=scan_calpha_traj(inputfile,AP.NATOM,N,M);
    else if (flag2=='d')
      ab=scandtraj(inputfile,numdihed,M);

    if (flag2=='d')
      FourierTransformdTraj(numdihed,M);
    else
      FourierTransformTraj(N, M);
    
    if (ns >= num_of_split_ini && ns < num_of_split_fin)
      for (i=0;i<M;++i)
	abF_ave[i]+=abF[i];
  }

  for (i=0;i<M;++i)
    abF_ave[i]=abF_ave[i]/(num_of_split_fin-num_of_split_ini);
  
  if (flag=='v') {
    ABS=0.0;
    for (i=0;i<M/2;++i) {
      ABS+=abF_ave[i];
    }
    if((outputfile=fopen(outputfilename,"w"))==NULL) {
      printf("There is not %s\n",outputfilename);
      exit(1);
    }
    for (i=0;i<M;++i)
      fprintf(outputfile,"%lf %8.3lf \n", 
	      (double)i/M/deltat/c,abF_ave[i]/ABS*DOF);
    fclose(outputfile);
  }
  else {
    ABS=0.0;
    pi=acos(-1.0);
    for (i=0;i<M/2;++i) {
      ABS+=(2.0*pi*i/M/deltat)*(2.0*pi*i/M/deltat)*abF_ave[i];
    }
    if((outputfile=fopen(outputfilename,"w"))==NULL) {
      printf("There is not %s\n",outputfilename);
      exit(1);
    }
    for (i=0;i<M;++i) {
      fprintf(outputfile,"%lf %13.8lf \n",
	      (double)i/M/deltat/c,
	      (2.0*pi*i/M/deltat)*(2.0*pi*i/M/deltat)*abF_ave[i]/ABS*DOF*2);
    }
    fclose(outputfile);
  }
  free(abF);

  return 0;
}

int scantraj(FILE *inputfile,int numatom, int N) {
  int i,j,f;
  double d,abs;
  FILE *outputfile;

  traj = malloc(sizeof(double)*N*numatom*3);

  abs=0.0;
  for (f=0;f<N;++f)  {
    for (i=0;i<numatom;++i)  {
      for (j=0;j<3;++j) {
	fscanf(inputfile,"%lf",&d);
	traj[f*numatom*3+i*3+j] = d*sqrt(AP.AMASS[i]);
	abs += traj[f*numatom*3+i*3+j]*traj[f*numatom*3+i*3+j];
      }
    }
  }
  
  return abs;
}

int scan_calpha_traj(FILE *inputfile,int numatom, int numres, int N) {
  int i,j,f;
  double d,abs;
  FILE *outputfile;

  traj = malloc(sizeof(double)*N*numres*3);

  abs=0.0;
  for (f=0;f<N;++f)  {
    for (i=0;i<numatom;++i)  {
      for (j=0;j<3;++j) {
	fscanf(inputfile,"%lf",&d);
	if (strncmp(AP.IGRAPH[i],"CA\0",3) == 0) {
	  traj[f*numres*3+i*3+j] = d*sqrt(AP.AMASS[i]);
	  abs += traj[f*numres*3+i*3+j]*traj[f*numres*3+i*3+j];
	}
      }
    }
  }

  return abs;
}

int scandtraj(FILE *inputfile, /*double *traj,*/ int numdihed, int N) {
  int i,j,f;
  double d,abs;
  FILE *outputfile;

  traj = malloc(sizeof(double)*N*numdihed);

  abs=0.0;
  for (f=0;f<N;++f)  {
    for (i=0;i<numdihed;++i)  {
      fscanf(inputfile,"%lf",&d);
      traj[f*numdihed+i] = d;
      abs += traj[f*numdihed+i]*traj[f*numdihed+i];
    }
  }

  return abs;
}

void FourierTransformTraj(int natom, int nstep){
  int i,j,k,na;
  double pi;
  fftw_complex *in,*out;
  fftw_plan p;
  FILE *outputfile;

  abF=(double *)malloc(sizeof(double)*nstep);
  for (i=0;i<nstep;++i){
    abF[i]=0.0;
  }

  for (na=0;na<natom;++na) {
    in  = fftw_malloc(sizeof(fftw_complex)*nstep);
    out = fftw_malloc(sizeof(fftw_complex)*nstep);
    p   = fftw_plan_dft_1d(nstep,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

    for (i=0;i<nstep;++i){
      in[i][0]=0.0;    in[i][1]=0.0;
    }

    for (i=0;i<nstep;++i)
      for (j=0;j<3;++j)
	in[i][0]+=traj[i*natom*3+na*3+j];
    
    fftw_execute(p);
    
    for (i=0;i<nstep;++i){
      abF[i] += out[i][0]*out[i][0]+out[i][1]*out[i][1];
    }
    
    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
  }
  
}

void FourierTransformdTraj(int numdihed, int nstep){
  int i,j,k,na;
  double pi;
  fftw_complex *in,*out;
  fftw_plan p;
  FILE *outputfile;

  if((abF=(double *)malloc(sizeof(double)*nstep))==NULL) {
    printf("error: cannot allocate abF\n");
    exit(1);
  }
  for (i=0;i<nstep;++i){
    abF[i]=0.0;
  }

  for (na=0;na<numdihed;++na) {
    in  = fftw_malloc(sizeof(fftw_complex)*nstep);
    out = fftw_malloc(sizeof(fftw_complex)*nstep);
    p   = fftw_plan_dft_1d(nstep,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
    for (i=0;i<nstep;++i){
      in[i][0]=0.0;    in[i][1]=0.0;
    }
    for (i=0;i<nstep;++i)
      in[i][0]=traj[i*numdihed+na];
    
    fftw_execute(p);
    
    for (i=0;i<nstep;++i){
      abF[i] += out[i][0]*out[i][0]+out[i][1]*out[i][1];
    }
    
    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
  }  
}
