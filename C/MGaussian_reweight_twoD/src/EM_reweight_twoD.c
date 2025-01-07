
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EM_reweight_twoD.h"
#include "Simpson_integ_twoD.h"
#include "Gaussian_twoD.h"
#include "EF.h"

void E_step_twoD_rw(int num_sim, int *n, double ***x,                     // # of simulation, # of snapshots, data,
		    double *minx, double *maxx, int num_Simpson,          // parameters for Simpson integration 1
		    double **k_umbrella, double **x0,                     // parameters for Simpson integration 2
		    int K, double **nyu_k,double ***Sigma_k,double *pi_k, // parameters of MG
		    double ***gammak_ij,  double **z_ik, double *z_i,     // responsibilities
		    double pi,int periodicflag, double periodicity) {
  int i,j,k,l;
  double S,S2;
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

  for (i=0;i<num_sim;++i) {
    z_i[i]=0.0;
    for (j=0;j<K;++j) {
      S=Simpson_integ_C_twoD_Gaussian(num_Simpson,minx,maxx,
  				      k_umbrella[i],x0[i],nyu_k[j],Sigma_k[j],pi,
  				      periodicflag,periodicity);
      z_ik[i][j]=pi_k[j]*S;
      z_i[i]+=z_ik[i][j];
    }
  }

}

void M_step_twoD_rw(int num_sim, int *n,int n_sum, double ***x,           // # of simulation, # of snapshots, data,
		    double *minx, double *maxx, int num_Simpson,          // parameters for Simpson integration 1
		    double **k_umbrella, double **x0,                     // parameters for Simpson integration 2
		    int K, double **nyu_k,double ***Sigma_k,double *pi_k, // parameters of MG
		    double ***gammak_ij,  double **z_ik, double *z_i,     // responsibilities
		    double pi, int periodicflag, double periodicity) {
  int i,j,k;
  double *Nk;
  double ***z_ik2,****z_ik3,z;
  double *S1,**S2,S;
  double *denominator;
  double v[2];

  S1=(double *)gcemalloc(sizeof(double)*2);

  S2=(double **)gcemalloc(sizeof(double *)*2);
  for (i=0;i<2;++i) S2[i]=(double *)gcemalloc(sizeof(double)*2);

  z_ik2=(double ***)gcemalloc(sizeof(double **)*2);
  for (i=0;i<2;++i) {
    z_ik2[i]=(double **)gcemalloc(sizeof(double *)*num_sim);
    for (j=0;j<num_sim;++j) {
      z_ik2[i][j]=(double *)gcemalloc(sizeof(double)*K);
    }
  }
  z_ik3=(double ****)gcemalloc(sizeof(double ***)*2);
  for (i=0;i<2;++i) {
    z_ik3[i]=(double ***)gcemalloc(sizeof(double **)*2);
    for (j=0;j<2;++j) {
      z_ik3[i][j]=(double **)gcemalloc(sizeof(double *)*num_sim);
      for (k=0;k<num_sim;++k) {
	z_ik3[i][j][k]=(double *)gcemalloc(sizeof(double)*K);
      }
    }
  }
  Nk=(double *)gcemalloc(sizeof(double)*K);


  for (i=0;i<num_sim;++i) {
    for (j=0;j<K;++j) {
      S1=Simpson_integ_C_x_uk_twoD_Gaussian(num_Simpson,minx,maxx,
					    k_umbrella[i],x0[i],nyu_k[j],Sigma_k[j],
					    pi,periodicflag,periodicity);
      z_ik2[0][i][j]=pi_k[j]*S1[0];
      z_ik2[1][i][j]=pi_k[j]*S1[1];
    }
  }

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
	nyu_k[k][0]+=(gammak_ij[k][i][j]*x[0][i][j]-z_ik2[0][i][k]/z_i[i])/(Nk[k]);
	nyu_k[k][1]+=(gammak_ij[k][i][j]*x[1][i][j]-z_ik2[1][i][k]/z_i[i])/(Nk[k]);
      }
    }
  }

  for (i=0;i<num_sim;++i) {
    for (j=0;j<K;++j) {
      S2=Simpson_integ_C_x_uk2_twoD_Gaussian(num_Simpson,minx,maxx,
					     k_umbrella[i],x0[i],nyu_k[j],Sigma_k[j],
					     pi,periodicflag,periodicity);
      z_ik3[0][0][i][j]=pi_k[j]*S2[0][0];
      z_ik3[0][1][i][j]=pi_k[j]*S2[0][1];
      z_ik3[1][0][i][j]=pi_k[j]*S2[1][0];
      z_ik3[1][1][i][j]=pi_k[j]*S2[1][1];
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

	Sigma_k[k][0][0]+=(gammak_ij[k][i][j]*v[0]*v[0]-z_ik3[0][0][i][k]/z_i[i])/(Nk[k]);
	Sigma_k[k][0][1]+=(gammak_ij[k][i][j]*v[0]*v[1]-z_ik3[0][1][i][k]/z_i[i])/(Nk[k]);
	Sigma_k[k][1][0]+=(gammak_ij[k][i][j]*v[1]*v[0]-z_ik3[1][0][i][k]/z_i[i])/(Nk[k]);
	Sigma_k[k][1][1]+=(gammak_ij[k][i][j]*v[1]*v[1]-z_ik3[1][1][i][k]/z_i[i])/(Nk[k]);
      }
    }
  }

  for (k=0;k<K;++k) {
    S=Simpson_integ_twoD_Gaussian(num_Simpson,minx,maxx,nyu_k[k],Sigma_k[k],pi);
    pi_k[k]=1.0/S/K;
  }
}

double EM_L_twoD_rw(int num_sim, int *n, double ***x,             // # of simulation, # of snapshots, data,
		    double *z_i,                                  // free energy differences bet. simulations
		    int K, double **nyu_k,double ***Sigma_k,double *pi_k, // parameters of MG
		    double pi) {
  int i,j,k;
  double L=0,lnp,X[2];

  for (i=0;i<num_sim;++i) {
    for (j=0;j<n[i];++j) {
      lnp=0.0;
      X[0]=x[0][i][j];
      X[1]=x[1][i][j];
      for (k=0;k<K;++k) {
	lnp+=1.0/z_i[i]*pi_k[k]*twoD_Gaussian(X,nyu_k[k],Sigma_k[k],pi);
      }
      lnp=log(lnp);
      L+=lnp;
    }
  }

  return L;
}
