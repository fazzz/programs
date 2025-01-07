
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "RAND.h"
#include "BOXMULL.h"
#include "MD.h"
#include "stringOTF.h"
#include "TOPO.h"
#include "PTL.h"

double String_cTheta_FASYS(double *theta,double *crd) {
  int i,j,k,t;
  double atom[2][4][3];
  double pi;

  pi=acos(-1.0);

  for (i=0;i<3;++i) {
    atom[0][0][i]=crd[0*3+i];
    atom[0][1][i]=crd[1*3+i];
    atom[0][2][i]=crd[2*3+i];
    atom[0][3][i]=crd[3*3+i];

    atom[1][0][i]=crd[1*3+i];
    atom[1][1][i]=crd[2*3+i];
    atom[1][2][i]=crd[3*3+i];
    atom[1][3][i]=crd[4*3+i];
  }

  theta[0]=dih(atom[0][0],atom[0][1],atom[0][2],atom[0][3]);
  theta[1]=dih(atom[1][0],atom[1][1],atom[1][2],atom[1][3]);

}

double String_Propagetor_Iso_FASYS(double *crd,double *vel,double *mass,int numatom,double IsoCoff,double dt,double *KE,double *PE,struct potential e,struct force f,double *z,double kappa,int numcv) {
  int i,j,k;
  double *velp,*crd_h,*frc,*crddx;
  double dx=1.0e-10;
  double sum;
  double dth_dx,*theta,*thetadx;

  velp=(double *)gcemalloc(sizeof(double)*numatom*3);
  crd_h=(double *)gcemalloc(sizeof(double)*numatom*3);
  crddx=(double *)gcemalloc(sizeof(double)*numatom*3);
  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  theta=(double *)gcemalloc(sizeof(double)*numcv);
  thetadx=(double *)gcemalloc(sizeof(double)*numcv);

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd_h[i*3+j]=crd[i*3+j]+dt*vel[i*3+j];

  ffL_calcffandforce(crd_h,numatom,&e,&f);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      frc[i*3+j]=-f.f_b[i*3+j]+f.f_a[i*3+j]+f.f_d[i*3+j];
      for (k=0;k<numatom*3;++k) crddx[k]=crd[k];
      crddx[i*3+j]+=dx;
      String_cTheta_FASYS(theta,crd);
      String_cTheta_FASYS(thetadx,crddx);
      for (k=0;k<numcv;++k) frc[i*3+j]-=kappa*(theta[k]-z[k])*(thetadx[k]-theta[k])/dx;
    }
  }
  *(PE)=e.p_b_t+e.p_a_t+e.p_d_t;

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) velp[i*3+j]=vel[i*3+j]+0.5*dt/mass[i]*frc[i*3+j];

  sum=0.0;
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) sum+=velp[i*3+j]*velp[i*3+j];
  sum=sqrt(sum);

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      vel[i*3+j]=IsoCoff*velp[i*3+j]/sum;
      crd[i*3+j]+=dt*vel[i*3+j];
    }
  }

  *(KE)=0.0;
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) *(KE)+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
}

double String_Propagetor(double *z_p,double *z,int numcv,double **M,double *theta,double kappa,double gamma,double dt){
  int i,j;

  for (i=0;i<numcv;++i) {
    z_p[i]=z[i];
    for (j=0;j<numcv;++j) 
      z_p[i]+=-kappa/gamma*dt*M[i][j]*(z[j]-theta[j]);
  }
}

double String_Repara(double **z,double **z_p,int numreplica,int numcv){
  int i,j,k;
  int *q;
  double *L,*l;
  double f;

  L=(double *)gcemalloc(sizeof(double)*numcv);
  l=(double *)gcemalloc(sizeof(double)*numcv);
  q=(int *)gcemalloc(sizeof(int)*numcv);

  L[0]=0.0;
  f=0.0;
  for (i=1;i<numreplica;++i) {
    for (j=0;j<numcv;++j) f+=(z_p[i][j]-z_p[i-1][j])*(z_p[i][j]-z_p[i-1][j]);
    L[i]=sqrt(f);
  }

  for (i=1;i<numreplica;++i) l[i]=L[numreplica-1]*i/(numreplica-1);

  for (i=1;i<numreplica;++i) {
    for (j=1;j<numreplica;++j) {
      if (L[j-1] < l[i] && l[i] <= L[j]) {
	q[i]=j;
	break;
      }      
    }
  }

  for (j=0;j<numcv;++j) {
      z[0][j]=z_p[0][j];
      z[numreplica-1][j]=z_p[numreplica-1][j];
  }

  for (i=1;i<numreplica-1;++i) {
    f=0.0;
    for (j=0;j<numcv;++j) f+=z_p[q[i]][j]-z_p[q[i-1]][j]*z_p[q[i]][j]-z_p[q[i-1]][j];
    f=sqrt(f);
    for (j=0;j<numcv;++j) z[i][j]=z_p[q[i-1]][j]+(l[i]-L[q[i-1]]-1.0)*(z_p[q[i]][j]-z_p[q[i-1]][j])/f;
  }
}


double String_cM_FASYS(double *crd,double *mass,double **M, int numatom,int numcv){
  int i,j,k,l,m;
  double dx=1.0e-10;
  double dth_dx;
  double *crddx;
  double *theta,*thetadx;
  double dxi,dxj;
  double pi;

  pi=acos(-1.0);
  theta=(double *)gcemalloc(sizeof(double)*2);
  thetadx=(double *)gcemalloc(sizeof(double)*2);
  crddx=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numcv;++i) {
    for (j=0;j<numcv;++j) { 
      M[i][j]=0.0;
      for (m=0;m<numatom*3;++m) crddx[m]=crd[m];
      for (k=0;k<numatom;++k) {
	for (l=0;l<3;++l) {
	  crddx[k*3+l]+=dx;
	  
	  String_cTheta_FASYS(theta,crd);
	  String_cTheta_FASYS(thetadx,crddx);
	  dxi=thetadx[i]-theta[i];
	  dxj=thetadx[i]-theta[i];
	  if (dxi>pi) dxi-=2.0*pi;
	  if (dxi<-1.0*pi) dxi+=2.0*pi;
	  if (dxj>pi) dxj-=2.0*pi;
	  if (dxj<-1.0*pi) dxj+=2.0*pi;
	  M[i][j]+=(dxi)/dx*(dxj)/dx/mass[k];
	}
      }
    }
  }
}

