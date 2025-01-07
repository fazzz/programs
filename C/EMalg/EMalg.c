
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EMalg.h"
#include "Gaussian.h"
#include "EF.h"

int E_step(int N,int K,double **x_n,double **nyu_k,double ***Sigma_k,double *pi_k, double **gamma_nk,double pi) {
  int i,j,k;
  double numerator,dinomirator;

  for (i=0;i<N;++i) {
    for (j=0;j<K;++j) {
      numerator=pi_k[j]*twoD_Gaussian(x_n[i],nyu_k[j],Sigma_k[j],pi);
      dinomirator=0.0;
      for (k=0;k<K;++k) {
	dinomirator+=pi_k[k]*twoD_Gaussian(x_n[i],nyu_k[k],Sigma_k[k],pi);
      }
      gamma_nk[i][j]=numerator/dinomirator;
    }
  }

}

double M_step(int N,int K,double **x_n,double **nyu_k,double ***Sigma_k,double *pi_k, double **gamma_nk) {
  int i,j,k;
  double *N_k;
  double vec[2];

  N_k=(double *)gcemalloc(sizeof(double)*K);

  for (i=0;i<K;++i) {
    N_k[i]=0.0;
    for (j=0;j<N;++j) {
      N_k[i]+=gamma_nk[j][i];
    }
  }

  for (i=0;i<K;++i) {
    pi_k[i]=N_k[i]/N;
  }

  for (i=0;i<K;++i) {
    nyu_k[i][0]=0.0;
    nyu_k[i][1]=0.0;
    for (j=0;j<N;++j) {
      nyu_k[i][0]+=gamma_nk[j][i]*x_n[j][0]/N_k[i];
      nyu_k[i][1]+=gamma_nk[j][i]*x_n[j][1]/N_k[i];
    }
  }

  for (i=0;i<K;++i) {

    Sigma_k[i][0][0]=0.0;
    Sigma_k[i][0][1]=0.0;
    Sigma_k[i][1][0]=0.0;
    Sigma_k[i][1][1]=0.0;

    for (j=0;j<N;++j) {
      vec[0]=x_n[j][0]-nyu_k[i][0];
      vec[1]=x_n[j][1]-nyu_k[i][1];

      Sigma_k[i][0][0]+=gamma_nk[j][i]*vec[0]*vec[0]/N_k[i];
      Sigma_k[i][0][1]+=gamma_nk[j][i]*vec[0]*vec[1]/N_k[i];
      Sigma_k[i][1][0]+=gamma_nk[j][i]*vec[1]*vec[0]/N_k[i];
      Sigma_k[i][1][1]+=gamma_nk[j][i]*vec[1]*vec[1]/N_k[i];
    }
    /*********************************************************/
    /* if (Sigma_k[i][0][0]==0.0 || Sigma_k[i][1][1]==0.0) { */
    /*   printf("i=%3d\n",i);				     */
    /* }						     */
    /*********************************************************/
  }

}

double EM_lnrou(int N,int K,double **x_n,double **nyu_k,double ***Sigma_k,double *pi_k, double pi) {
  int i,j,k;
  double lnrou=0.0,sum=0.0;

  for (i=0;i<N;++i) {
    sum=0.0;
    for (j=0;j<K;++j) {
      sum+=pi_k[j]*twoD_Gaussian(x_n[i],nyu_k[j],Sigma_k[j],pi);
    }
    lnrou+=log(sum);
  }

  return lnrou;
}

int E_step_fprob(int Nx, double minx,double dx,int Ny,double miny,double dy,
		 int K,double **prob,double **nyu_k,double ***Sigma_k,double *pi_k, double ***gamma_nk,
		 double pi) {
  int i,j,k,n=0,nx,ny;
  double *x;
  double numerator,dinomirator;

  x=(double *)gcemalloc(sizeof(double)*2);

  n=0;
  for (nx=0;nx<Nx;++nx) {
    x[0]=minx+nx*dx;
    for (ny=0;ny<Ny;++ny) {
      x[1]=miny+ny*dy;

      for (i=0;i<K;++i) {
	/********************************/
        /* if (i==16) {		        */
	/*   printf("i=%d\n",i);        */
	/* }			        */
        /********************************/
	numerator=pi_k[i]*twoD_Gaussian(x,nyu_k[i],Sigma_k[i],pi);
	dinomirator=0.0;
	for (j=0;j<K;++j) {
	  dinomirator+=pi_k[j]*twoD_Gaussian(x,nyu_k[j],Sigma_k[j],pi);
	}
	/***************************************/
        /* if (i==16) {			       */
	/*   printf("%lf\n",numerator);	       */
	/* }				       */
        /***************************************/
	gamma_nk[nx][ny][i]=numerator/dinomirator;
      }
    }
  }

}

double M_step_fprob(int Nx,double minx,double dx,int Ny, double miny,double dy,
		    /*int N,*/int K,double **prob,double **nyu_k,double ***Sigma_k,double *pi_k, double ***gamma_nk) {
  int i,j,k,n=0,nx,ny;
  double *x;
  double *N_k;
  double vec[2];

  x=(double *)gcemalloc(sizeof(double)*2);

  N_k=(double *)gcemalloc(sizeof(double)*K);

  for (i=0;i<K;++i) {
    N_k[i]=0.0;
    for (nx=0;nx<Nx;++nx) {
      for (ny=0;ny<Ny;++ny) {
	if (prob[nx][ny]!=0.0) // debug
	  N_k[i]+=gamma_nk[nx][ny][i]*prob[nx][ny];
      }
    }
  }

  for (i=0;i<K;++i) {
    pi_k[i]=N_k[i]/*/N*/;
  }

  for (i=0;i<K;++i) {
    nyu_k[i][0]=0.0;
    nyu_k[i][1]=0.0;
    n=0;
    for (nx=0;nx<Nx;++nx) {
      x[0]=minx+nx*dx;
      for (ny=0;ny<Ny;++ny) {
	if (prob[nx][ny]!=0.0) { // debug
	  x[1]=miny+ny*dy;

	  /**************************************************************/
	  /* if (i==16 /\*&& gamma_nk[nx][ny][i]!=0*\/ ) {	      */
	  /*   printf("gamma_nk=( %lf )\n",gamma_nk[nx][ny][i]);	      */
	  /* }							      */
	  /**************************************************************/

	  nyu_k[i][0]+=gamma_nk[nx][ny][i]*x[0]/N_k[i]*prob[nx][ny];
	  nyu_k[i][1]+=gamma_nk[nx][ny][i]*x[1]/N_k[i]*prob[nx][ny];
	}
      }
    }
  }

  for (i=0;i<K;++i) {

    Sigma_k[i][0][0]=0.0;
    Sigma_k[i][0][1]=0.0;
    Sigma_k[i][1][0]=0.0;
    Sigma_k[i][1][1]=0.0;


    n=0;
    for (nx=0;nx<Nx;++nx) {
      x[0]=minx+nx*dx;
      for (ny=0;ny<Ny;++ny) {
	if (prob[nx][ny]!=0.0) { // debug
	  x[1]=miny+ny*dy;

	  vec[0]=x[0]-nyu_k[i][0];
	  vec[1]=x[1]-nyu_k[i][1];

	  /*****************************************************************************************/
	  /* if (i==16 && gamma_nk[nx][ny][i]!=0.0)						 */
	  /*   printf("gamma_nk=%lf, vec=( %lf %lf )\n",gamma_nk[nx][ny][i],vec[0],vec[1]);	 */
	  /*****************************************************************************************/

	  Sigma_k[i][0][0]+=gamma_nk[nx][ny][i]*vec[0]*vec[0]/N_k[i]*prob[nx][ny];
	  Sigma_k[i][0][1]+=gamma_nk[nx][ny][i]*vec[0]*vec[1]/N_k[i]*prob[nx][ny];
	  Sigma_k[i][1][0]+=gamma_nk[nx][ny][i]*vec[1]*vec[0]/N_k[i]*prob[nx][ny];
	  Sigma_k[i][1][1]+=gamma_nk[nx][ny][i]*vec[1]*vec[1]/N_k[i]*prob[nx][ny];
	}
      }
    }
  }

  for (i=0;i<K;++i) {                // debug
    if (Sigma_k[i][0][0] <= 1.0e-6 )
      Sigma_k[i][0][0]=0.1;
    if (Sigma_k[i][1][1] <= 1.0e-6 )
      Sigma_k[i][1][1]=0.1;
  }

}

double EM_lnrou_fprob(int Nx, double minx,double dx, int Ny,double miny,double dy,
		      /*int N,*/int K,double **prob,double **nyu_k,double ***Sigma_k,double *pi_k, double pi) {
  int i,j,k,n=0,nx,ny;
  double *x;
  double tG;
  double lnrou=0.0,sum=0.0;

  x=(double *)gcemalloc(sizeof(double)*2);

  for (nx=0;nx<Nx;++nx) {
    x[0]=minx+nx*dx;
    for (ny=0;ny<Ny;++ny) {
      if (prob[nx][ny]!=0.0) { // debug
	x[1]=miny+ny*dy;
	sum=0.0;
	for (j=0;j<K;++j) {
	  tG=twoD_Gaussian(x,nyu_k[j],Sigma_k[j],pi)*prob[nx][ny];
	  sum+=pi_k[j]*tG;
	}
	if (sum!=0.0)
	  lnrou+=log(sum);
      }
    }
  }

  return lnrou;
}

double *Create_mixed_twoD_GaussianMap(double minx,double maxx,double dx,
				      double miny,double maxy,double dy,
				      double **nyu,double ***Sigma, double *pi_k, int K, double pi){
  int i,j,k;
  int n=0;
  double *x;
  double *prob;
  
  x=(double *)gcemalloc(sizeof(double)*2);
  prob=(double *)gcemalloc(sizeof(double)*1);

  x[0]=minx;
  for (i=0;x[0]<maxx;++i){
    x[0]=minx+dx*i;
    x[1]=miny;
    for (j=0;x[1]<maxy;++j){
      x[1]=miny+dy*j;
      prob[n]=0.0;
      for (k=0;k<K;++k){
	prob[n]+=pi_k[k]*twoD_Gaussian(x,nyu[k],Sigma[k],pi);
      }
      ++n;
      prob=(double *)gcerealloc(prob,sizeof(double)*n);
    }
  }

  return prob;
}

