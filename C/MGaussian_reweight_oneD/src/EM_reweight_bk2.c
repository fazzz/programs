
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
      //      S2=Simpson_integ_oneD_Gaussian(num_Simpson,minx,maxx,pi);
      z_ik[i][j]=pi_k[j]*S;
      z_i[i]+=z_ik[i][j];
    }
  }

  /********************************************************************************************************/
  /* for (i=0;i<num_sim;++i) {										  */
  /*   for (j=0;j<n[i];++j) {										  */
  /*     delta=fabs(x[i][j]-x0[i]);									  */
  /*     if (periodicflag==ON) if (fabs(delta)>0.5*periodicity) delta=/\*2.0**\/periodicity-delta;	  */
  /*     //      w[i][j]=z_i[i]*exp(0.5*k_umbrella[i]*delta*delta/k_B/T);				  */
  /*     w[i][j]=/\*1.0/\*\/z_i[i]*exp(-0.5*k_umbrella[i]*delta*delta/k_B/T);				  */
  /*     /\************************************************\/						  */
  /*     /\* if (delta>pi) 				      *\/					  */
  /*     /\* 	printf("%d %d %8.6lf\n",i,j,delta);   *\/						  */
  /*     /\* printf("%d %d %8.6lf\n",i,j,w[i][j]);	      *\/					  */
  /*     /\************************************************\/						  */
  /*     for (k=0;k<K;++k) {										  */
  /* 	gammak_ij[k][i][j]=w[i][j]*gammak_ij[k][i][j];							  */
  /*     }												  */
  /*   }												  */
  /* }													  */
  /********************************************************************************************************/

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

  /*********************************************************/
  /* denominator=(double *)gcemalloc(sizeof(double)*K);	   */
  z_ik2=(double **)gcemalloc(sizeof(double *)*num_sim);
  z_ik3=(double **)gcemalloc(sizeof(double *)*num_sim);
  for (i=0;i<num_sim;++i) {
    z_ik2[i]=(double *)gcemalloc(sizeof(double)*K);
    z_ik3[i]=(double *)gcemalloc(sizeof(double)*K);
  }
  Nk1=(double *)gcemalloc(sizeof(double)*K);
  Nk2=(double *)gcemalloc(sizeof(double)*K);		   
  /*********************************************************/

  /***************************************************************/
  /* z=0.0;							 */
  /* for (i=0;i<K;++i) {					 */
  /*   S=Simpson_integ_oneD_Gaussian(num_Simpson,minx,maxx, pi); */
  /*   z+=pi_k[i]*S;						 */
  /* }								 */
  /***************************************************************/

  for (i=0;i<num_sim;++i) {
    for (j=0;j<K;++j) {
      /**************************************************************************************/
      /* S1=Simpson_integ_C_x_oneD_Gaussian(num_Simpson,minx,maxx,			    */
      /* 					 k_umbrella[i],x0[i],nyu_k[j],Sigma_k[j],   */
      /* 					 pi,periodicflag,periodicity);		    */
      /**************************************************************************************/
      S1=Simpson_integ_C_x_uk_oneD_Gaussian(num_Simpson,minx,maxx,
					    k_umbrella[i],x0[i],nyu_k[j],Sigma_k[j],
					    pi,periodicflag,periodicity);
      z_ik2[i][j]=pi_k[j]*S1;
    }
  }

  for (k=0;k<K;++k) {
    Nk1[k]=0.0;
    //    Nk2[k]=0.0;
    for (i=0;i<num_sim;++i) {
      for (j=0;j<n[i];++j) {
	Nk1[k]+=gammak_ij[k][i][j];
	//	Nk2[k]+=z_ik[i][k]/z_i[i];
      }
    }
  }

  /***********************/
  /* sum=0.0;		 */
  /* for (i=0;i<K;++i) { */
  /*   sum+=Nk1[i];	 */
  /* }			 */
  /* 			 */
  /* sum=0.0;		 */
  /* for (i=0;i<K;++i) { */
  /*   sum+=Nk2[i];	 */
  /* }			 */
  /***********************/

  for (k=0;k<K;++k) {
    nyu_k[k]=0.0;
    for (i=0;i<num_sim;++i) {
      for (j=0;j<n[i];++j) {
	//	nyu_k[k]+=(gammak_ij[k][i][j]*x[i][j]/*+z_ik2[i][k]/z_i[i]*/)/(Nk1[k]/*+Nk2[k]*/);
	//	nyu_k[k]+=(gammak_ij[k][i][j]*x[i][j]-z_ik2[i][k]/z_i[i])/(Nk1[k]-Nk2[k]);
	//	nyu_k[k]+=(gammak_ij[k][i][j]*x[i][j]-z_ik2[i][k]/z_i[i])/(Nk1[k]-Nk2[k]);
	//	nyu_k[k]+=(gammak_ij[k][i][j]*x[i][j])/(Nk1[k]);
	nyu_k[k]+=(gammak_ij[k][i][j]*x[i][j]-z_ik2[i][k]/z_i[i])/(Nk1[k]);
      }
    }
  }

  for (i=0;i<num_sim;++i) {
    for (j=0;j<K;++j) {
      /***************************************************************************************/
      /* S2=Simpson_integ_C_x2_oneD_Gaussian(num_Simpson,minx,maxx,			     */
      /* 					  k_umbrella[i],x0[i],nyu_k[j],Sigma_k[j],   */
      /* 					  pi,periodicflag,periodicity);		     */
      /***************************************************************************************/
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
	v=x[i][j]-nyu_k[k];
	//	if (periodicflag==ON) if (fabs(v)>0.5*periodicity) v=2.0*periodicity-v;
	//	Sigma_k[k]+=(gammak_ij[k][i][j]*v*v/*+z_ik3[i][k]/z_i[i]*/)/(Nk1[k]/*+Nk2[k]*/);
	//	Sigma_k[k]+=(gammak_ij[k][i][j]*v*v-z_ik3[i][k]/z_i[i])/(Nk1[k]-Nk2[k]);
	//	Sigma_k[k]+=(gammak_ij[k][i][j]*v*v-z_ik3[i][k]/z_i[i])/(Nk1[k]-Nk2[k]);
	//	Sigma_k[k]+=(gammak_ij[k][i][j]*v*v)/(Nk1[k]);
	Sigma_k[k]+=(gammak_ij[k][i][j]*v*v-z_ik3[i][k]/z_i[i])/(Nk1[k]);
      }
    }
  }

  for (k=0;k<K;++k) {
    Sigma_k[k]=sqrt(Sigma_k[k]);
  }

  /******************************************************************************************/
  /* for (k=0;k<K;++k) {								    */
  /*   denominator[k]=0.0;								    */
  /*   for (i=0;i<num_sim;++i) {							    */
  /*     S=Simpson_integ_C_oneD_Gaussian(num_Simpson,minx,maxx,				    */
  /* 				      k_umbrella[i],x0[i],nyu_k[k],Sigma_k[k],pi,	    */
  /* 				      periodicflag,periodicity);			    */
  /*     denominator[k]+=n[i]*S/z_i[i];							    */
  /*   }										    */
  /* }											    */
  /* 											    */
  /* sum=0.0;										    */
  /* for (k=0;k<K;++k) {								    */
  /*   sum+=denominator[k];								    */
  /* }											    */
  /******************************************************************************************/

  /********************************************************************************************************/
  /* sum=0.0;												  */
  /* for (i=0;i<num_sim;++i) {										  */
  /*   for (j=0;j<n[i];++j) {										  */
  /*     delta=fabs(x[i][j]-x0[i]);									  */
  /*     if (periodicflag==ON) if (fabs(delta)>0.5*periodicity) delta=/\*2.0**\/periodicity-delta;	  */
  /*     for (k=0;k<K;++k) {										  */
  /* 	//  	sum+=z_i[i]*exp(0.5*k_umbrella[i]*delta*delta/k_B/T);					  */
  /* 	sum+=1.0/z_i[i]*exp(-0.5*k_umbrella[i]*delta*delta/k_B/T);					  */
  /*     }												  */
  /*   }												  */
  /* }													  */
  /********************************************************************************************************/

  /*****************************************/
  /* sum=0.0;				   */
  /* for (i=0;i<num_sim;++i) {		   */
  /*   for (j=0;j<n[i];++j) {		   */
  /*     for (k=0;k<K;++k) {		   */
  /* 	sum+=gammak_ij[k][i][j];	   */
  /*     }				   */
  /*   }				   */
  /* }					   */
  /*****************************************/

  for (k=0;k<K;++k) {
    S=Simpson_integ_oneD_Gaussian(num_Simpson,minx,maxx,nyu_k[k],Sigma_k[k],pi);
    pi_k[k]=1.0/S/K;
  }

  for (k=0;k<K;++k) {
    //    pi_k[k]=(Nk1[k]/*+Nk2[k]*/)/n_sum/*/2.0*/;
    //    pi_k[k]=(Nk1[k]+Nk2[k])/n_sum/2.0;
    //    pi_k[k]=Nk1[k]/n_sum;
    //    pi_k[k]=Nk1[k]/sum;
    /**********************************/
    /* sum=0.0;			      */
    /* for (i=0;i<num_sim;++i) {      */
    /*   sum+=n[i]/z_i[i]*z_ik[i][k]; */
    /* }			      */
    /* pi_k[k]=Nk1[k]/sum;	      */
    /**********************************/
    //    pi_k[k]=Nk1[k]/denominator[k];
    //    pi_k[k]=Nk1[k]/Nk2[k]*pi_k[k];
    ;
  }

  /***********************/
  /* sum=0.0;		 */
  /* for (k=0;k<K;++k) { */
  /*   sum+=pi_k[k];	 */
  /* }			 */
  /***********************/
  
  /* for (k=0;k<K;++k) {    */
  /*   pi_k[k]=pi_k[k]/sum; */
  /* }			    */
  /**************************/

  /******************************************************************************************/
  /* for (i=0;i<num_sim;++i) {								    */
  /*   for (j=0;j<K;++j) {								    */
  /*     S1=Simpson_integ_C_x_oneD_Gaussian(num_Simpson,minx,maxx,			    */
  /*     					 k_umbrella[i],x0[i],nyu_k[j],Sigma_k[j],   */
  /*     					 pi,periodicflag,periodicity);		    */
  /*     z_ik2[i][j]=pi_k[j]*S1;							    */
  /*   }										    */
  /* }											    */
  /******************************************************************************************/

  /********************************************/
  /* for (k=0;k<K;++k) {		      */
  /*   Nk1[k]=0.0;			      */
  /*   Nk2[k]=0.0;			      */
  /*   for (i=0;i<num_sim;++i) {	      */
  /*     for (j=0;j<n[i];++j) {		      */
  /* 	Nk1[k]+=gammak_ij[k][i][j];	      */
  /* 	Nk2[k]+=z_ik[i][k]/z_i[i];	      */
  /*     }				      */
  /*   }				      */
  /* }					      */
  /********************************************/

  /************************************/
  /* for (k=0;k<K;++k) {	      */
  /*   pi_k[k]=Nk1[k]/Nk2[k]*pi_k[k]; */
  /* }				      */
  /************************************/

  /**************************************************/
  /* for (i=0;i<K;++i) {			    */
  /*   printf("%d: %8.3e %8.3e\n",i,Nk1[i],Nk2[i]); */
  /* }						    */
  /* printf("\n");				    */
  /**************************************************/

}

double EM_L_rw(int num_sim, int *n, double **x,                    // # of simulation, # of snapshots, data,
	       double *z_i,                                        // free energy differences bet. simulations
	       int K, double *nyu_k,double *Sigma_k,double *pi_k,  // parameters of MG
	       double pi) {
  int i,j,k;
  double L,L1=0.0,L2=0.0,lnp;

  /****************************************************************************/
  /* for (i=0;i<num_sim;++i) {						      */
  /*   for (j=0;j<n[i];++j) {						      */
  /*     lnp=0.0;							      */
  /*     for (k=0;k<K;++k) {						      */
  /* 	lnp+=pi_k[k]*oneD_Gaussian(x[i][j],nyu_k[k],Sigma_k[k],pi);	      */
  /*     }								      */
  /*     lnp=log(lnp);							      */
  /*     L1+=lnp;							      */
  /*   }								      */
  /* }									      */
  /* 									      */
  /* for (i=0;i<num_sim;++i) {						      */
  /*   L2+=n[i]*log(z_i[i]);						      */
  /* }									      */
  /****************************************************************************/

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

  //  L=L1;
  //  L=L1-L2;
  //  L=L1-L2;

  //  printf("%5.3e %5.3e\n",L1,L2);

  return L;
}
