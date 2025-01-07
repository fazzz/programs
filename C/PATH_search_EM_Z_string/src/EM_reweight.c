
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EM_reweight.h"
#include "Gaussian.h"
#include "Simpson_integ.h"
#include "EF.h"

void E_step_oneD(int num_sim, int *n, double **x,                      // # of simulation, # of snapshots, data,
		 double minx, double maxx, int num_Simpson,            // parameters for Simpson integration 1
		 double *k_umbrella, double *x0,                       // parameters for Simpson integration 2
		 int K, double *nyu_k,double *Sigma_k,double *pi_k,    // parameters of MG
		 double ***gammak_ij,  double **inv_f_ik, double *f_i, // responsibilities
		 double pi) {
  int i,j,k,l;
  double S;
  double numerator,dinomirator;

  for (i=0;i<num_sim;++i) {
    for (j=0;j<n[i];++j) {
      for (k=0;k<K;++k) {
	numerator=pi_k[k]*oneD_Gaussian(x[i][j],nyu_k[k],Sigma_k[k],pi);
	dinomirator=0.0;
	for (l=0;l<K;++l) 
	  dinomirator+=pi_k[l]*oneD_Gaussian(x[i][j],nyu_k[l],Sigma_k[l],pi);
	gammak_ij[k][i][j]=numerator/dinomirator;
      }
    }
  }

  for (i=0;i<num_sim;++i) {
    f_i[i]=0.0;
    for (j=0;j<K;++j) {
      S=Simpson_integ_C_oneD_Gaussian(num_Simpson,minx,maxx,k_umbrella[j],x0[j],nyu_k[j],Sigma_k[j],pi);
      inv_f_ik[i][j]=pi_k[j]*S;
      f_i[i]+=inv_f_ik[i][j];
    }
    f_i[i]=1.0/f_i[i];
  }
}

void M_step_oneD(int num_sim, int *n,int n_sum, double **x,            // # of simulation, # of snapshots, data,
		 double minx, double maxx, int num_Simpson,            // parameters for Simpson integration 1
		 double *k_umbrella, double *x0,                       // parameters for Simpson integration 2
		 int K, double *nyu_k,double *Sigma_k,double *pi_k,    // parameters of MG
		 double ***gammak_ij,  double **inv_f_ik, double *f_i, // responsibilities
		 double pi) {
  int i,j,k;
  double *Nk1,*Nk2;
  double **f_ik2,**f_ik3;
  double S1,S2;
  double dinomirator1,dinomirator2;

  f_ik2=(double **)gcemalloc(sizeof(double *)*num_sim);
  f_ik3=(double **)gcemalloc(sizeof(double *)*num_sim);
  for (i=0;i<num_sim;++i) {
    f_ik2[i]=(double *)gcemalloc(sizeof(double)*K);
    f_ik3[i]=(double *)gcemalloc(sizeof(double)*K);
  }
  Nk1=(double *)gcemalloc(sizeof(double)*K);
  Nk2=(double *)gcemalloc(sizeof(double)*K);

  for (i=0;i<num_sim;++i) {
    for (j=0;j<K;++j) {
      S1=Simpson_integ_C_x_oneD_Gaussian(num_Simpson,minx,maxx,k_umbrella[j],x0[j],nyu_k[j],Sigma_k[j],pi);
      S2=Simpson_integ_C_x2_oneD_Gaussian(num_Simpson,minx,maxx,k_umbrella[j],x0[j],nyu_k[j],Sigma_k[j],pi);
      f_ik2[i][j]=pi_k[j]*S1;
      f_ik3[i][j]=pi_k[j]*S2;
    }
  }

  for (i=0;i<K;++i) {
    Nk1[i]=0.0; Nk2[i]=0.0;
  }

  for (i=0;i<num_sim;++i)
    for (j=0;j<n[i];++j) 
      for (k=0;k<K;++k) 
	Nk1[k]+=gammak_ij[k][i][j];

  for (i=0;i<num_sim;++i)
    for (j=0;j<K;++j)
      Nk2[j]+=n[i]*inv_f_ik[i][j]*f_i[i];

  for (i=0;i<K;++i) {
    nyu_k[i]=0.0;
    Sigma_k[i]=0.0;
    pi_k[i]=(Nk1[i]+Nk2[i])/n_sum/2.0;
  }

  for (k=0;k<K;++k) {
    for (i=0;i<num_sim;++i) {
      for (j=0;j<n[i];++j) {
	nyu_k[k]+=(gammak_ij[k][i][j]*x[i][j]+f_ik2[i][k]*f_i[i])/(Nk1[k]+Nk2[k]);

	Sigma_k[k]+=(gammak_ij[k][i][j]*(x[i][j]-nyu_k[k])*(x[i][j]-nyu_k[k])+f_ik3[i][k]*f_i[i])/(Nk1[k]+Nk2[k]);
      }
    }
  }
}

double EM_L(int num_sim, int *n, double **x,                    // # of simulation, # of snapshots, data,
	    double *f,                                          // free energy differences bet. simulations
	    int K, double *nyu_k,double *Sigma_k,double *pi_k,  // parameters of MG
	    double pi) {
  int i,j,k;
  double L,L1=0.0,L2=0.0,lnp;

  for (i=0;i<num_sim;++i) {
    for (j=0;j<n[i];++j) {
      lnp=0.0;
      for (k=0;k<K;++k) {
	lnp+=(pi_k[k]*oneD_Gaussian(x[i][j],nyu_k[k],Sigma_k[k],pi));
      }
      lnp=log(lnp);
      L1+=lnp;
    }
  }

  for (i=0;i<num_sim;++i) {
    L2+=n[i]*log(f[i]);
  }

  L=L1-L2;

  return L;
}
