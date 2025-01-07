
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EM.h"
//#include "Gaussian.h"
#include "Simpson_integ.h"
#include "EF.h"

void E_step_oneD(int num_sim, int *n, double **x,                      // # of simulation, # of snapshots, data,
		 double minx, double maxx, int num_Simpson,            // parameters for Simpson integration 1
		 int K, double *nyu_k,double *Sigma_k,double *pi_k,    // parameters of MG
		 double ***gammak_ij,                                  // responsibilities
		 double pi,int periodicflag, double periodicity) {
  int i,j,k,l;
  double S;
  double numerator,dinomirator;

  double f;

  for (i=0;i<num_sim;++i) {
    for (j=0;j<n[i];++j) {
      for (k=0;k<K;++k) {
	numerator=pi_k[k]*oneD_Gaussian(x[i][j],nyu_k[k],Sigma_k[k],pi);
	dinomirator=0.0; for (l=0;l<K;++l) dinomirator+=pi_k[l]*oneD_Gaussian(x[i][j],nyu_k[l],Sigma_k[l],pi);
	gammak_ij[k][i][j]=numerator/dinomirator;
      }
    }
  }

}

void M_step_oneD(int num_sim, int *n,int n_sum, double **x,            // # of simulation, # of snapshots, data,
		 double minx, double maxx, int num_Simpson,            // parameters for Simpson integration 1
		 int K, double *nyu_k,double *Sigma_k,double *pi_k,    // parameters of MG
		 double ***gammak_ij,                                  // responsibilities
		 double pi, int periodicflag, double periodicity) {
  int i,j,k;
  double *Nk1,*Nk2;
  double **z_ik2,**z_ik3;
  double S1,S2;
  double v;

  Nk1=(double *)gcemalloc(sizeof(double)*K);

  for (k=0;k<K;++k) {
    Nk1[k]=0.0;
    for (i=0;i<num_sim;++i) {
      for (j=0;j<n[i];++j) {
	Nk1[k]+=gammak_ij[k][i][j];
      }
    }
  }

  for (k=0;k<K;++k) {
    nyu_k[k]=0.0;
    for (i=0;i<num_sim;++i) {
      for (j=0;j<n[i];++j) {
	nyu_k[k]+=(gammak_ij[k][i][j]*x[i][j])/(Nk1[k]);
      }
    }
  }

  for (k=0;k<K;++k) {
    Sigma_k[k]=0.0;
    for (i=0;i<num_sim;++i) {
      for (j=0;j<n[i];++j) {
	v=x[i][j]-nyu_k[k];
	if (periodicflag==ON) if (fabs(v)>0.5*periodicity) v=2.0*periodicity-v;
	Sigma_k[k]+=(gammak_ij[k][i][j]*v*v)/(Nk1[k]);
      }
    }
  }

  for (k=0;k<K;++k) {
    Sigma_k[k]=sqrt(Sigma_k[k]);
  }

  for (k=0;k<K;++k) {
    pi_k[k]=(Nk1[k])/n_sum;
  }

}

double EM_L(int num_sim, int *n, double **x,                    // # of simulation, # of snapshots, data,
	    int K, double *nyu_k,double *Sigma_k,double *pi_k,  // parameters of MG
	    double pi) {
  int i,j,k;
  double L,L1=0.0,L2=0.0,lnp;

  for (i=0;i<num_sim;++i) {
    for (j=0;j<n[i];++j) {
      lnp=0.0;
      for (k=0;k<K;++k) {
	lnp+=pi_k[k]*oneD_Gaussian(x[i][j],nyu_k[k],Sigma_k[k],pi);
      }
      lnp=log(lnp);
      L1+=lnp;
    }
  }

  L=L1;

  return L;
}
