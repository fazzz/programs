
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EM_reweight.h"
//#include "Gaussian.h"
#include "Simpson_integ.h"
#include "EF.h"

void E_step_oneD(int num_sim, int *n, double **x,                      // # of simulation, # of snapshots, data,
		 double minx, double maxx, int num_Simpson,            // parameters for Simpson integration 1
		 double *k_umbrella, double *x0,                       // parameters for Simpson integration 2
		 int K, double *nyu_k,double *Sigma_k,double *pi_k,    // parameters of MG
		 double ***gammak_ij,  double **z_ik, double *z_i,     // responsibilities
		 double pi,int periodicflag, double periodicity) {
  int i,j,k,l;
  double S;
  double numerator,dinomirator;

  double f;

  for (i=0;i<num_sim;++i) {
    for (j=0;j<n[i];++j) {
      for (k=0;k<K;++k) {
	numerator=pi_k[k]*oneD_Gaussian(x[i][j],nyu_k[k],Sigma_k[k],pi);
	/*************************************************************/
        /* if (numerator==0.0) {				     */
	/*   f=oneD_Gaussian(x[i][j],nyu_k[k],Sigma_k[k],pi);	     */
	/*   printf("debug");					     */
	/* }							     */
        /*************************************************************/
	dinomirator=0.0; for (l=0;l<K;++l) dinomirator+=pi_k[l]*oneD_Gaussian(x[i][j],nyu_k[l],Sigma_k[l],pi);
	gammak_ij[k][i][j]=numerator/dinomirator;
      }
    }
  }

  for (i=0;i<num_sim;++i) {
    z_i[i]=0.0;
    for (j=0;j<K;++j) {
      S=Simpson_integ_C_oneD_Gaussian(num_Simpson,minx,maxx,
				      k_umbrella[/*j*/i],x0[/*j*/i],nyu_k[j],Sigma_k[j],pi,
				      periodicflag,periodicity);
      z_ik[i][j]=pi_k[j]*S;
      z_i[i]+=z_ik[i][j];
    }
  }
}

void M_step_oneD(int num_sim, int *n,int n_sum, double **x,            // # of simulation, # of snapshots, data,
		 double minx, double maxx, int num_Simpson,            // parameters for Simpson integration 1
		 double *k_umbrella, double *x0,                       // parameters for Simpson integration 2
		 int K, double *nyu_k,double *Sigma_k,double *pi_k,    // parameters of MG
		 double ***gammak_ij,  double **z_ik, double *z_i,     // responsibilities
		 double pi, int periodicflag, double periodicity) {
  int i,j,k;
  double *Nk1,*Nk2;
  double **z_ik2,**z_ik3;
  double S1,S2;
  double dinomirator,deltax_nyu;
  double v;

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
      S1=Simpson_integ_C_x_oneD_Gaussian(num_Simpson,minx,maxx,
					 k_umbrella[/*j*/i],x0[/*j*/i],nyu_k[j],Sigma_k[j],
					 pi,periodicflag,periodicity);
      S2=Simpson_integ_C_x2_oneD_Gaussian(num_Simpson,minx,maxx,
					  k_umbrella[/*j*/i],x0[/*j*/i],nyu_k[j],Sigma_k[j],
					  pi,periodicflag,periodicity);
      z_ik2[i][j]=pi_k[j]*S1;
      z_ik3[i][j]=pi_k[j]*S2;
    }
  }

  /*****************************/
  /* for (i=0;i<K;++i) {       */
  /*   Nk1[i]=0.0; Nk2[i]=0.0; */
  /* }			       */
  /*****************************/

  for (k=0;k<K;++k) {
    Nk1[k]=0.0;
    Nk2[k]=0.0;
    for (i=0;i<num_sim;++i) {
      for (j=0;j<n[i];++j) {
	Nk1[k]+=gammak_ij[k][i][j];
	Nk2[k]+=z_ik[i][/*j*/k]/z_i[i];
      }
    }
  }

  /***********************/
  /* for (i=0;i<K;++i) { */
  /*   nyu_k[i]=0.0;	 */
  /*   Sigma_k[i]=0.0;	 */
  /* }			 */
  /***********************/

  for (k=0;k<K;++k) {
    nyu_k[k]=0.0;
    for (i=0;i<num_sim;++i) {
      for (j=0;j<n[i];++j) {
	//	nyu_k[k]+=(gammak_ij[k][i][j]*x[i][j]+z_ik2[i][k]/z_i[i])/(Nk1[k]+Nk2[k]);
	nyu_k[k]+=(gammak_ij[k][i][j]*x[i][j])/Nk1[k];

	//	deltax_nyu=fabs(x[i][j]-nyu_k[k]);
	//	if (periodicflag==ON) if (fabs(deltax_nyu)>0.5*periodicity) deltax_nyu=2.0*periodicity-deltax_nyu;
	//	Sigma_k[k]+=(gammak_ij[k][i][j]*deltax_nyu*deltax_nyu+z_ik3[i][k]/z_i[i])/(Nk1[k]+Nk2[k]);
	//	Sigma_k[k]+=(gammak_ij[k][i][j]*(x[i][j]-nyu_k[k])*(x[i][j]-nyu_k[k]))/(Nk1[k]);
	/*********************************/
        /* if (Sigma_k[k]==0.0) {	 */
	/*   printf("debug");		 */
	/* }				 */
        /*********************************/
      }
    }
    //    pi_k[i]=(Nk1[i]+Nk2[i])/n_sum/2.0;
    //    pi_k[k]=(Nk1[k])/n_sum;
  }

  /********************************************************************************************************/
  /* for (k=0;k<K;++k) {										  */
  /*   nyu_k[k]=nyu_k[k]/Nk1[k];									  */
  /* }													  */
  /* 													  */
  for (k=0;k<K;++k) {
    Sigma_k[k]=0.0;
    for (i=0;i<num_sim;++i) {
      for (j=0;j<n[i];++j) {
	v=x[i][j]-nyu_k[k];
  	Sigma_k[k]+=(gammak_ij[k][i][j]/**v*v*/)/Nk1[k];
      }
    }
  }
  /* 													  */
  /* for (k=0;k<K;++k) {										  */
  /*   Sigma_k[k]=Sigma_k[k]/Nk1[k];									  */
  /* }													  */
  /********************************************************************************************************/

  for (k=0;k<K;++k) {
    pi_k[k]=(Nk1[k])/n_sum;
  }

}

double EM_L(int num_sim, int *n, double **x,                    // # of simulation, # of snapshots, data,
	    double *z_i,                                        // free energy differences bet. simulations
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

  for (i=0;i<num_sim;++i) {
    L2+=n[i]*log(z_i[i]);
  }

  //  L=L1+L2;
  L=/*-*/L1;

  return L;
}
