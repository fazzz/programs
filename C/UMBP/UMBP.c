
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#include "UMBP.h"
#include "PTL.h"

#include "TOPO.h"
#include "LA.h"
#include "mymath.h"

double UMB_calc_dihetype_ff(double *crd,int numatom,int *pairp,int num,double *fcp,double *dih_equ,double *p,double **f){
  int i,j,k,l;
  int ii,jj,kk,ll;

  double fi[3],fj[3],fk[3],fl[3];
  double m[3],n[3],m_n[3],n_n[3],lm,ln;
  double vij[3],vkj[3],vkl[3];
  double lkj;
  double vijvkj,vklvkj;

  double judge;

  double dvdpsi;
  
  double atom[4][3];
  double dihed,delta;
  double p_t=0.0;

  double pi;

  pi=acos(-1.0);

  for (i=0;i<num;++i) p[i] = 0.0;
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) f[i][j] = 0.0;
  
  for (i=0;i<num;++i) {
    ii=pairp[i*4+0];
    jj=pairp[i*4+1];
    kk=pairp[i*4+2];
    ll=pairp[i*4+3];

    for (j=0;j<3;++j) {
      atom[0][j]=crd[ii*3+j];
      atom[1][j]=crd[jj*3+j];
      atom[2][j]=crd[kk*3+j];
      atom[3][j]=crd[ll*3+j];
    }

    for (j=0;j<3;++j) {
      vij[j] = atom[1][j]-atom[0][j];
      vkj[j] = atom[1][j]-atom[2][j];
      vkl[j] = atom[3][j]-atom[2][j];
    }
    lkj=sqrt(inprod(vkj,vkj,3));
    
    outprod(vij,vkj,m);
    outprod(vkj,vkl,n);
    lm=sqrt(inprod(m,m,3));
    ln=sqrt(inprod(n,n,3));
    for (j=0;j<3;++j) {
      m_n[j]=m[j]/lm;
      n_n[j]=n[j]/ln;
    }
    
    dihed=inprod(m_n,n_n,3);
    if (dihed>=1.0)
      dihed=0.0;
    else if (dihed<=-1.0)
      dihed=pi;
    else
      dihed=acos(dihed);
    if (inprod(vij,n,3)>0) dihed=-dihed;
    if (dihed<0.0) dihed=2.0*pi+dihed;
    
    vijvkj=inprod(vij,vkj,3);
    vklvkj=inprod(vkl,vkj,3);
    
    if ((delta=dihed-dih_equ[i])>pi) delta-=2.0*pi;
    else if ((delta=dihed-dih_equ[i])<-1.0*pi) delta+=2.0*pi;

    dvdpsi=2.0*fcp[i]*delta;
    
    p[i] = fcp[i]*(delta)*(delta);
    p_t += p[i];
    
    for (j=0;j<3;++j) {
      fi[j] = -dvdpsi*lkj*m[j]/(lm*lm);
      fl[j] =  dvdpsi*lkj*n[j]/(ln*ln);
      fj[j] = (-fi[j]+(vijvkj/(lkj*lkj))*fi[j]-(vklvkj/(lkj*lkj))*fl[j]);
      fk[j] = (-fl[j]-(vijvkj/(lkj*lkj))*fi[j]+(vklvkj/(lkj*lkj))*fl[j]);
      
      f[ii][j] += fi[j]*4.184070*100.0;
      f[jj][j] += fj[j]*4.184070*100.0;
      f[kk][j] += fk[j]*4.184070*100.0;
      f[ll][j] += fl[j]*4.184070*100.0;
    }
  }

  return p_t;
}

