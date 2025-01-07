
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "EMalg_non_square_range.h"
#include "Gaussian.h"
#include "EF.h"

int E_step_fprob_nsr(int Nx, double minx,double dx,int Ny,double miny,double dy,
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
      if (prob[nx][ny]!=0.0) {

	for (i=0;i<K;++i) {
	  numerator=pi_k[i]*twoD_Gaussian(x,nyu_k[i],Sigma_k[i],pi);
	  dinomirator=0.0;
	  for (j=0;j<K;++j) {
	    dinomirator+=pi_k[j]*twoD_Gaussian(x,nyu_k[j],Sigma_k[j],pi);
	  }
	  gamma_nk[nx][ny][i]=numerator/dinomirator;
	}
      }
    }
  }

}

double M_step_fprob_nsr(int Nx,double minx,double dx,int Ny, double miny,double dy,
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

double EM_lnrou_fprob_nsr(int Nx, double minx,double dx, int Ny,double miny,double dy,
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
	  //	  if (tG!=0.0) // 2013-06-07
	  sum+=pi_k[j]*tG;
	}
	if (sum!=0.0)
	  lnrou+=log(sum);
      }
    }
  }

  return lnrou;
}

double EM_lnrou_fprob_nsr2(int Nx, double minx,double dx, int Ny,double miny,double dy,
			   /*int N,*/int K,double **prob,double **nyu_k,double ***Sigma_k,double *pi_k, double pi) {
  int i,j,k,n=0,nx,ny;
  double *x;
  double tG,a;
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
	  if (tG!=0.0) {
	    a=pi_k[j]*tG;
	    sum+=a;
	  }
	}
	if (sum!=0.0)
	  lnrou+=log(sum);
      }
    }
  }

  return lnrou;
}

double M_step_fprob_nsr_2013_06_07(int Nx,double minx,double dx,int Ny, double miny,double dy,
				   /*int N,*/int K,double **prob,double **nyu_k,double ***Sigma_k,double *pi_k, double ***gamma_nk) {
  int i,j,k,n=0,nx,ny;
  double *x;
  double *N_k,a;
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

  for (i=0;i<K;++i) {
    if ((a=Sigma_k[i][0][0]*Sigma_k[i][1][1]-Sigma_k[i][1][0]*Sigma_k[i][0][1]) <= 1.0e-6) {
      Sigma_k[i][0][0]=Sigma_k[i][0][0]*10.0;
      Sigma_k[i][1][1]=Sigma_k[i][1][1]*10.0;
    }
  }

}
