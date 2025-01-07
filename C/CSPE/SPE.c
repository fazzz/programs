
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fftw3.h"

#include "SPE.h"

void FT(int natom, int nstep,double *trj, double *spe){
  int i,j,k,l;
  fftw_complex *in,*out;
  fftw_plan p;
  FILE *outputfile;

  double sumspe,K;

  for (i=0;i<nstep;++i)
    spe[i]=0.0;

  for (i=0;i<natom;++i) {
    for (j=0;j<3;++j) {
      in  = fftw_malloc(sizeof(fftw_complex)*nstep);
      out = fftw_malloc(sizeof(fftw_complex)*nstep);
      p   = fftw_plan_dft_1d(nstep,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

      for (k=0;k<nstep;++k) {
	in[k][0]+=trj[k*natom*3+i*3+j];
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
	in[j][0]+=dtrj[j*numdihed+i];
    
    fftw_execute(p);
    
    for (j=0;j<nstep;++j)
      spe[j] += out[j][0]*out[j][0]+out[j][1]*out[j][1];
    
    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
  }  
}

void FT_ts(int nstep,double *dat, double *spe){
  int i,j,k,l;
  fftw_complex *in,*out;
  fftw_plan p;

  double sumspe,K;

  for (i=0;i<nstep;++i) spe[i]=0.0;

  in  = fftw_malloc(sizeof(fftw_complex)*nstep);
  out = fftw_malloc(sizeof(fftw_complex)*nstep);
  p   = fftw_plan_dft_1d(nstep,in,out,FFTW_FORWARD,FFTW_ESTIMATE);

  for (i=0;i<nstep;++i) {
    in[i][0]=dat[i];
    in[i][1]=0.0;
    out[i][0]=0.0;
    out[i][1]=0.0;
  }
  fftw_execute(p);
    
  for (i=0;i<nstep;++i) spe[i] = out[i][0]*out[i][0]+out[i][1]*out[i][1];

  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
}

void CSP_decom(int natom, int nstep,double *trj, double *spe){
  int i,j,k,l;
  fftw_complex *in,*out;
  fftw_plan p;
  FILE *outputfile;

  double sumspe,K;

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
	spe[k*natom+i] += out[k][0]/nstep*out[k][0]/nstep+out[k][1]/nstep*out[k][1]/nstep;

      fftw_destroy_plan(p);
      fftw_free(in); fftw_free(out);
    }  
  }    
}

void CSPd_decom(int numdihed, int nstep,double *dtrj, double *spe){
  int i,j,k;
  fftw_complex *in,*out;
  fftw_plan p;
  FILE *outputfile;

  for (i=0;i<nstep*numdihed;++i)
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
      spe[j*numdihed+i] = out[j][0]*out[j][0]+out[j][1]*out[j][1];
    
    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
  }  
}
