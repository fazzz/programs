
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "K_means_non_square_range.h"
#include "EF.h"

double Kmeans_Estep_nsr(int N, int K, double **x, double **nyu, double **gamma_nk, double **dist){
  int i,j;
  int index_min;
  double min,sum;

  for (i=0;i<N;++i) {
    for (j=0;j<K;++j) {
      gamma_nk[i][j]=0.0;
    }
  }

  for (i=0;i<N;++i) {
    min=dist[i][0];
    index_min=0;
    for (j=0;j<K;++j) {
      if (min>dist[i][j]) {
	min=dist[i][j];
	index_min=j;
      }
    }
    gamma_nk[i][index_min]=1.0;
  }

  sum=0.0;
  for (i=0;i<N;++i) {
    for (j=0;j<K;++j) {
      sum+=gamma_nk[i][j];
    }
  }

  return sum;
}

double Kmeans_Mstep_nsr(int N, int K, double **x, double **nyu, double **gamma_nk){
  int i,j,k;
  double sum_gamma_k,sum_gamma_x_k[2];

  for (j=0;j<K;++j) {
    sum_gamma_k=0.0;
    for (i=0;i<2;++i) {
      sum_gamma_x_k[i]=0.0;
    }
    for (k=0;k<N;++k) {
      sum_gamma_k+=gamma_nk[k][j];
      for (i=0;i<2;++i) {
	sum_gamma_x_k[i]+=gamma_nk[k][j]*x[k][i];
      }
    }

    for (i=0;i<2;++i) {
      nyu[j][i]=sum_gamma_x_k[i]/sum_gamma_k;
    }
  }

  return 0.0;
}

double Kmeans_J_nsr(int N, int K,double **x,double **nyu, double **gamma_nk, double **dist){
  int i,j;
  double J;

  for (i=0;i<N;++i) {
    for (j=0;j<K;++j) {
      dist[i][j]=(x[i][0]-nyu[j][0])*(x[i][0]-nyu[j][0])+(x[i][1]-nyu[j][1])*(x[i][1]-nyu[j][1]);
    }
  }

  J=0.0;

  for (i=0;i<N;++i) {
    for (j=0;j<K;++j) {
      J+=gamma_nk[i][j]*dist[i][j];
    }
  }

  return J;
}

double Kmeans_Estep_fprob_nsr(int Nx, double minx,double dx,int Ny,double miny,double dy,
			      int K, double **prob, double **nyu, double ***gamma_nk, double ***dist){
  int i,j,k,n=0,nx,ny;
  double *x;
  int index_min;
  double min,sum;

  x=(double *)gcemalloc(sizeof(double)*2);

  for (nx=0;nx<Nx;++nx) {
    for (ny=0;ny<Ny;++ny) {
      if (prob[nx][ny]!=0.0) {
	for (j=0;j<K;++j) {
	  gamma_nk[nx][ny][j]=0.0;
	}
      }
    }
  }

  for (nx=0;nx<Nx;++nx) {
    for (ny=0;ny<Ny;++ny) {
      if (prob[nx][ny]!=0.0) {
	min=dist[nx][ny][0];
	index_min=0;
	for (j=0;j<K;++j) {
	  if (min>dist[nx][ny][j]) {
	    min=dist[nx][ny][j];
	    index_min=j;
	  }
	}
	gamma_nk[nx][ny][index_min]=1.0;
      }
    }
  }

  sum=0.0;
  for (nx=0;nx<Nx;++nx) {
    for (ny=0;ny<Ny;++ny) {
      if (prob[nx][ny]!=0.0) {
	for (j=0;j<K;++j) {
	  sum+=gamma_nk[nx][ny][j];
	}
      }
    }
  }

  return sum;
}

double Kmeans_Mstep_fprob_nsr(int Nx, double minx,double dx,int Ny,double miny,double dy,
			      int K, double **prob, double **nyu, double ***gamma_nk){
  int i,j,k,nx,ny;
  double *x;
  double sum_gamma_k,sum_gamma_x_k[2];

  x=(double *)gcemalloc(sizeof(double)*2);

  for (j=0;j<K;++j) {
    sum_gamma_k=0.0;
    for (i=0;i<2;++i) {
      sum_gamma_x_k[i]=0.0;
    }
    for (nx=0;nx<Nx;++nx) {
      x[0]=minx+nx*dx;
      for (ny=0;ny<Ny;++ny) {
	x[1]=miny+ny*dy;

	if (prob[nx][ny]!=0.0) {
	  sum_gamma_k+=gamma_nk[nx][ny][j]*prob[nx][ny];
	  for (i=0;i<2;++i) {
	    sum_gamma_x_k[i]+=gamma_nk[nx][ny][j]*x[i]*prob[nx][ny];
	  }
	}
      }
    }

    for (i=0;i<2;++i) {
      if (sum_gamma_k != 0.0) // 2013-06-07
	nyu[j][i]=sum_gamma_x_k[i]/sum_gamma_k;
    }
  }

  return 0.0;
}

double Kmeans_J_fprob_nsr(int Nx, double minx,double dx,int Ny,double miny,double dy,
			  int K,double **prob,double **nyu, double ***gamma_nk, double ***dist){
  int i,j,nx,ny;
  double *x;
  double J;

  x=(double *)gcemalloc(sizeof(double)*2);

  for (nx=0;nx<Nx;++nx) {
    x[0]=minx+nx*dx;
    for (ny=0;ny<Ny;++ny) {
      x[1]=miny+ny*dy;
      if (prob[nx][ny]!=0.0) {
	for (j=0;j<K;++j) {
	  dist[nx][ny][j]=(x[0]-nyu[j][0])*(x[0]-nyu[j][0])+(x[1]-nyu[j][1])*(x[1]-nyu[j][1]);
	}
      }
    }
  }

  J=0.0;

  for (nx=0;nx<Nx;++nx) {
    for (ny=0;ny<Ny;++ny) {
      if (prob[nx][ny]!=0.0) {
	for (j=0;j<K;++j) {
	  J+=gamma_nk[nx][ny][j]*dist[nx][ny][j]*prob[nx][ny];
	}
      }
    }
  }

  return J;
}
