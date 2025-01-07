
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EM_reweight.h"
#include "Simpson_integ.h"
#include "EF.h"

void E_step_oneD_rw(int num_sim, int *n, double **x,                      // # of simulation, # of snapshots, data,
		    double minx, double maxx, int num_Simpson,            // parameters for Simpson integration 1
		    double *k_umbrella, double *x0,                       // parameters for Simpson integration 2
		    int K, double *nyu_k,double *Sigma_k,double *pi_k,    // parameters of MG
		    double ***gammak_ij,  double **z_ik, double *z_i,     // responsibilities
		    double pi,int periodicflag, double periodicity) {
  int i,j,k,l;
  double S,S2;
  double numerator,dinomirator,delta;

  double f;

  double k_B=1.98723e-3;
  double T=300;

  double **w;

  w=(double **)gcemalloc(sizeof(double *)*num_sim);
  for (i=0;i<num_sim;++i) {
    w[i]=(double *)gcemalloc(sizeof(double)*n[i]);
  }

  for (i=0;i<num_sim;++i) {
    for (j=0;j<n[i];++j) {
      for (k=0;k<K;++k) {
	numerator=pi_k[k]*oneD_Gaussian(x[i][j],nyu_k[k],Sigma_k[k],pi);
	dinomirator=0.0; for (l=0;l<K;++l) dinomirator+=pi_k[l]*oneD_Gaussian(x[i][j],nyu_k[l],Sigma_k[l],pi);
	gammak_ij[k][i][j]=numerator/dinomirator;
      }
    }
  }

  for (i=0;i<num_sim;++i) {
    z_i[i]=0.0;
    for (j=0;j<K;++j) {
      S=Simpson_integ_C_oneD_Gaussian(num_Simpson,minx,maxx,
  				      k_umbrella[i],x0[i],nyu_k[j],Sigma_k[j],pi,
  				      periodicflag,periodicity);
      z_ik[i][j]=pi_k[j]*S;
      z_i[i]+=z_ik[i][j];
    }
  }

}

void M_step_oneD_rw(int num_sim, int *n,int n_sum, double **x,            // # of simulation, # of snapshots, data,
		    double minx, double maxx, int num_Simpson,            // parameters for Simpson integration 1
		    double *k_umbrella, double *x0,                       // parameters for Simpson integration 2
		    int K, double *nyu_k,double *Sigma_k,double *pi_k,    // parameters of MG
		    double ***gammak_ij,  double **z_ik, double *z_i,     // responsibilities
		    double pi, int periodicflag, double periodicity) {
  int i,j,k;
  double *Nk1,*Nk2;
  double **z_ik2,**z_ik3,z;
  double S1,S2,S;
  double *denominator;
  double v;
  double delta,sum;

  double k_B=1.98723e-3;
  double T=300;

  double *nyu_k_new;
  double *Sigma_k_new;

  int l;
  double numerator,dinomirator;

  nyu_k_new=(double *)gcemalloc(sizeof(double)*num_sim);
  Sigma_k_new=(double *)gcemalloc(sizeof(double)*num_sim);

  z_ik2=(double **)gcemalloc(sizeof(double *)*num_sim);
  z_ik3=(double **)gcemalloc(sizeof(double *)*num_sim);
  for (i=0;i<num_sim;++i) {
    z_ik2[i]=(double *)gcemalloc(sizeof(double)*K);
    z_ik3[i]=(double *)gcemalloc(sizeof(double)*K);
  }
  Nk1=(double *)gcemalloc(sizeof(double)*K);
  Nk2=(double *)gcemalloc(sizeof(double)*K);		   

  for (i=0;i<num_sim;++i) {
    for (j=0;j<K;++j) {
      S1=Simpson_integ_C_x_uk_oneD_Gaussian(num_Simpson,minx,maxx,
					    k_umbrella[i],x0[i],nyu_k[j],Sigma_k[j],
					    pi,periodicflag,periodicity);
      z_ik2[i][j]=pi_k[j]*S1;
    }
  }

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
	nyu_k[k]+=(gammak_ij[k][i][j]*x[i][j]-z_ik2[i][k]/z_i[i])/(Nk1[k]);
      }
    }
  }

  /***************************************************/
  /* if (periodicflag==ON) {			     */
  /*   for (k=0;k<K;++k) {			     */
  /*     while (nyu_k[k]</\*-3.0*pi*\/minx) {	     */
  /* 	nyu_k[k]+=periodicity;			     */
  /*     }					     */
  /*     while (nyu_k[k]>/\*3.0*pi*\/maxx) {	     */
  /* 	nyu_k[k]-=periodicity;			     */
  /*     }					     */
  /*   }					     */
  /* }						     */
  /***************************************************/

  for (i=0;i<num_sim;++i) {
    for (j=0;j<K;++j) {
      S2=Simpson_integ_C_x_uk2_oneD_Gaussian(num_Simpson,minx,maxx,
					     k_umbrella[i],x0[i],nyu_k[j],Sigma_k[j],
					     pi,periodicflag,periodicity);
      z_ik3[i][j]=pi_k[j]*S2;
    }
  }

  for (k=0;k<K;++k) {
    Sigma_k[k]=0.0;
    for (i=0;i<num_sim;++i) {
      for (j=0;j<n[i];++j) {
	v=fabs(x[i][j]-nyu_k[k]);
	//	if (periodicflag==ON) if (v>0.5*periodicity) v=/*2.0**/periodicity-v;
	Sigma_k[k]+=(gammak_ij[k][i][j]*v*v-z_ik3[i][k]/z_i[i])/(Nk1[k]);
      }
    }
  }

  for (k=0;k<K;++k) {
    Sigma_k[k]=sqrt(Sigma_k[k]);
  }

  /**********************************************************************************/
  /* for (k=0;k<K;++k) {							    */
  /*   S=Simpson_integ_oneD_Gaussian(num_Simpson,minx,maxx,nyu_k[k],Sigma_k[k],pi); */
  /*   pi_k[k]=1.0/S/K;								    */
  /* }										    */
  /**********************************************************************************/

  /***********************************************************************************/
  /* S=0.0;									     */
  /* for (k=0;k<K;++k) {							     */
  /*   S+=Simpson_integ_oneD_Gaussian(num_Simpson,minx,maxx,nyu_k[k],Sigma_k[k],pi); */
  /* }										     */
  /* for (k=0;k<K;++k) {							     */
  /*   pi_k[k]=1.0/S;								     */
  /* }										     */
  /***********************************************************************************/
  /***********************************************************************************/
  /* S=0.0;									     */
  /* for (k=0;k<K;++k) {							     */
  /*   S+=Simpson_integ_oneD_Gaussian(num_Simpson,minx,maxx,nyu_k[k],Sigma_k[k],pi); */
  /* }										     */
  /* for (k=0;k<(int)(K/2);++k) {						     */
  /*   pi_k[k]=0.9/S;								     */
  /* }										     */
  /* for (k=(int)(K/2);k<K;++k) {						     */
  /*   pi_k[k]=0.1/S;								     */
  /* }										     */
  /***********************************************************************************/
  S=0.0;
  for (k=0;k<(int)(K/4);++k) {
    S+=Simpson_integ_oneD_Gaussian(num_Simpson,minx,maxx,nyu_k[k],Sigma_k[k],pi);
  }
  for (k=0;k<(int)(K/4);++k) {
    pi_k[k]=0.60/S;
  }
  S=0.0;
  for (k=(int)(K/4);k<(int)(K/4*2);++k) {
    S+=Simpson_integ_oneD_Gaussian(num_Simpson,minx,maxx,nyu_k[k],Sigma_k[k],pi);
  }
  for (k=(int)(K/4);k<(int)(K/4*2);++k) {
    pi_k[k]=0.30/S;
  }
  S=0.0;
  for (k=(int)(K/4*2);k<(int)(K/4*3);++k) {
    S+=Simpson_integ_oneD_Gaussian(num_Simpson,minx,maxx,nyu_k[k],Sigma_k[k],pi);
  }
  for (k=(int)(K/4*2);k<(int)(K/4*3);++k) {
    pi_k[k]=0.15/S;
  }
  S=0.0;
  for (k=(int)(K/4*3);k<K;++k) {
    S+=Simpson_integ_oneD_Gaussian(num_Simpson,minx,maxx,nyu_k[k],Sigma_k[k],pi);
  }
  for (k=(int)(K/4*3);k<K;++k) {
    pi_k[k]=0.15/S;
  }

}

double EM_L_rw(int num_sim, int *n, double **x,                    // # of simulation, # of snapshots, data,
	       double *z_i,                                        // free energy differences bet. simulations
	       int K, double *nyu_k,double *Sigma_k,double *pi_k,  // parameters of MG
	       double pi) {
  int i,j,k;
  double L=0.0,L1=0.0,L2=0.0,lnp;

  for (i=0;i<num_sim;++i) {
    for (j=0;j<n[i];++j) {
      lnp=0.0;
      for (k=0;k<K;++k) {
	lnp+=1.0/z_i[i]*pi_k[k]*oneD_Gaussian(x[i][j],nyu_k[k],Sigma_k[k],pi);
      }
      lnp=log(lnp);
      L+=lnp;
    }
  }

  return L;
}
