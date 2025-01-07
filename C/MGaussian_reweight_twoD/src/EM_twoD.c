
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EM_twoD.h"
#include "Gaussian_twoD.h"
#include "EF.h"

void E_step_twoD(int num_sim, int *n, double ***x,                     // # of simulation, # of snapshots, data,
		 int K, double **nyu_k,double ***Sigma_k,double *pi_k, // parameters of MG
		 double ***gammak_ij,                                  // responsibilities
		 double pi) {
  int i,j,k,l;
  double numerator,dinomirator;
  double X[2];

  for (i=0;i<num_sim;++i) {
    for (j=0;j<n[i];++j) {
      for (k=0;k<K;++k) {
	X[0]=x[0][i][j];
	X[1]=x[1][i][j];

	numerator=pi_k[k]*twoD_Gaussian(X,nyu_k[k],Sigma_k[k],pi);
	dinomirator=0.0; for (l=0;l<K;++l) dinomirator+=pi_k[l]*twoD_Gaussian(X,nyu_k[l],Sigma_k[l],pi);
	gammak_ij[k][i][j]=numerator/dinomirator;
      }
    }
  }
}

void M_step_twoD(int num_sim, int *n,int n_sum, double ***x,           // # of simulation, # of snapshots, data,
		 int K, double **nyu_k,double ***Sigma_k,double *pi_k, // parameters of MG
		 double ***gammak_ij,                                  // responsibilities
		 double pi) {
  int i,j,k;
  double *Nk;
  double v[2];

  Nk=(double *)gcemalloc(sizeof(double)*K);

  for (k=0;k<K;++k) {
    Nk[k]=0.0;
    for (i=0;i<num_sim;++i) {
      for (j=0;j<n[i];++j) {
	Nk[k]+=gammak_ij[k][i][j];
      }
    }
  }

  for (k=0;k<K;++k) {
    nyu_k[k][0]=0.0;
    nyu_k[k][1]=0.0;
    for (i=0;i<num_sim;++i) {
      for (j=0;j<n[i];++j) {
	nyu_k[k][0]+=(gammak_ij[k][i][j]*x[0][i][j])/Nk[k];
	nyu_k[k][1]+=(gammak_ij[k][i][j]*x[1][i][j])/Nk[k];
      }
    }
  }

  for (k=0;k<K;++k) {
    Sigma_k[k][0][0]=0.0;
    Sigma_k[k][0][1]=0.0;
    Sigma_k[k][1][0]=0.0;
    Sigma_k[k][1][1]=0.0;

    for (i=0;i<num_sim;++i) {
      for (j=0;j<n[i];++j) {
	v[0]=x[0][i][j]-nyu_k[k][0];
	v[1]=x[1][i][j]-nyu_k[k][1];

	Sigma_k[k][0][0]+=(gammak_ij[k][i][j]*v[0]*v[0])/Nk[k];
	Sigma_k[k][0][1]+=(gammak_ij[k][i][j]*v[0]*v[1])/Nk[k];
	Sigma_k[k][1][0]+=(gammak_ij[k][i][j]*v[1]*v[0])/Nk[k];
	Sigma_k[k][1][1]+=(gammak_ij[k][i][j]*v[1]*v[1])/Nk[k];
      }
    }
  }

  for (k=0;k<K;++k) pi_k[k]=Nk[k]/n_sum;

}

double EM_L_twoD(int num_sim, int *n, double ***x,                      // # of simulation, # of snapshots, data,
		 int K, double **nyu_k,double ***Sigma_k,double *pi_k,  // parameters of MG
		 double pi) {
  int i,j,k;
  double X[2];
  double L=0.0,lnp;

  for (i=0;i<num_sim;++i) {
    for (j=0;j<n[i];++j) {
      X[0]=x[0][i][j];
      X[1]=x[1][i][j];
      lnp=0.0;
      for (k=0;k<K;++k) {
	lnp+=pi_k[k]*twoD_Gaussian(X,nyu_k[k],Sigma_k[k],pi);
      }
      lnp=log(lnp);
      L+=lnp;
    }
  }

  return L;
}
