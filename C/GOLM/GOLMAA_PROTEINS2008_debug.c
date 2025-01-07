
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLMAA_PROTEINS2008_debug.h"

#include "TOPO.h"
#include "LA.h"
#include "mymath.h"

#include "PTL.h"

#define ON 1
#define OFF 0

double GOLMAA_PROTEINS2008_debug_ff_calcff(double *crd, int numatom,struct potential_GOLMAA_PROTEINS2008 *ene, int bflag, int aflag, int dflag, int nflag ) {
  int i,j;

  if (bflag==ON)
    (*ene).p_b_t=GOLMAA_PROTEINS2008_debug_ff_calcBOND(crd,numatom,(*ene).p_b,(*ene).f_b,(*ene).Kb,(*ene).bon_equ,(*ene).pairs_bond,(*ene).num_bond);

  if (aflag==ON)
    (*ene).p_a_t=GOLMAA_PROTEINS2008_debug_ff_calcANGLE(crd,numatom,(*ene).p_a,(*ene).f_a,(*ene).Ka,(*ene).ang_equ,(*ene).pairs_angl,(*ene).num_angl);

  if (dflag==ON)
    (*ene).p_d_t=GOLMAA_PROTEINS2008_debug_ff_calcDIHE(crd,numatom,(*ene).p_d,(*ene).f_d,(*ene).Kd1,(*ene).Kd2,(*ene).dih_equ,(*ene).pairs_dihe,(*ene).num_dihe);

  if (nflag==ON)
    GOLMAA_PROTEINS2008_debug_ff_calcff_nonlocal(crd,numatom,
				     (*ene).NC_index,(*ene).numNC,(*ene).NotNC_index,(*ene).numNotNC,
				     (*ene).ALJ_natatt,(*ene).BLJ_natatt,(*ene).ep_natatt,(*ene).ALJ_repul,
				     (*ene).p_natatt,(*ene).p_repul,(*ene).f_natatt,(*ene).f_repul);

  (*ene).p_t=0.0;
  (*ene).p_natatt_t=0.0;
  (*ene).p_repul_t=0.0;
  if (nflag==ON) {
    for (i=0;i<(*ene).numNC;++i) (*ene).p_natatt_t+=(*ene).p_natatt[i];
    for (i=0;i<(*ene).numNotNC;++i) (*ene).p_repul_t+=(*ene).p_repul[i];
  }

  if (bflag==ON)
    (*ene).p_t+=(*ene).p_b_t;
  if (aflag==ON)
    (*ene).p_t+=(*ene).p_a_t;
  if (dflag==ON)
    (*ene).p_t+=(*ene).p_d_t;

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j)
      (*ene).f_t[i][j]=0.0;

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      if (bflag==ON)
	(*ene).f_t[i][j]+=(*ene).f_b[i][j];
      if (aflag==ON)
	(*ene).f_t[i][j]+=(*ene).f_a[i][j];
      if (dflag==ON)
	(*ene).f_t[i][j]+=(*ene).f_d[i][j];
      if (nflag==ON)
	(*ene).f_t[i][j]+=(*ene).f_natatt[i][j]+(*ene).f_repul[i][j];
    }
  }

  return (*ene).p_t;
}

double GOLMAA_PROTEINS2008_debug_ff_calcBOND(double *crd,int numatom,double *p_b,double **f_b,double Kb,double *bon_equ,int **pairs,int numbond){
  int i,j,k;
  int ii,jj;
  double f;
  double lenij;
  double atom[2][3];
  double p_b_t=0.0;

  for (i=0;i<numbond;++i) p_b[i] = 0.0;
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) f_b[i][j] = 0.0;
  
  for (i=0;i<numbond;++i) {
    ii=pairs[i][0];
    jj=pairs[i][1];
    for (j=0;j<3;++j) {
      atom[0][j]=crd[ii*3+j];
      atom[1][j]=crd[jj*3+j];
    }
  
    lenij = len(atom[0],atom[1]);
    p_b[i]=Kb*(lenij-bon_equ[i])*(lenij-bon_equ[i]);
    p_b_t += p_b[i];
    for (j=0;j<3;++j) {
      f = 2.0*Kb*(lenij-bon_equ[i])*(atom[1][j]-atom[0][j])/lenij*4.184070*100.0;
      f_b[ii][j] += f;
      f_b[jj][j] += -f;
    }
  }

  return p_b_t;
}

double GOLMAA_PROTEINS2008_debug_ff_calcANGLE(double *crd,int numatom,double *p_a,double **f_a,double Ka,double *ang_equ,int **pairs,int numangl){
  int i,j,k,l;
  int ii,jj,kk;
  double atom[3][3];
  double lenij,lenkj;
  double vij[3],vkj[3];
  double cosijk,angijk;
  double f1,f2;
  double p_a_t=0.0;

  for (i=0;i<numangl;++i) p_a[i] = 0.0;
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) f_a[i][j] = 0.0;

  for (i=0;i<numangl;++i) {
    ii=pairs[i][0];
    jj=pairs[i][1];
    kk=pairs[i][2];
    for (j=0;j<3;++j) {
      atom[0][j]=crd[ii*3+j];
      atom[1][j]=crd[jj*3+j];
      atom[2][j]=crd[kk*3+j];
    }

    lenij = len(atom[0],atom[1]);
    lenkj = len(atom[2],atom[1]);
    for (j=0;j<3;++j) {
      vij[j]=atom[1][j]-atom[0][j];
      vkj[j]=atom[1][j]-atom[2][j];
    }
    cosijk=inprod(vij,vkj,3);
    cosijk=cosijk/lenij/lenkj;
    angijk = acos(cosijk);

    angijk = ang(atom[0],atom[1],atom[2]);

    p_a[i] = Ka*(angijk-ang_equ[i])*(angijk-ang_equ[i]);
    p_a_t += p_a[i];
    
    for (j=0;j<3;++j) {
      f1 = -2.0*Ka*(angijk-ang_equ[i])/(lenij*sin(angijk))*(vkj[j]/lenkj-cosijk*vij[j]/lenij)*4.184070*100.0;
      f2 = -2.0*Ka*(angijk-ang_equ[i])/(lenkj*sin(angijk))*(vij[j]/lenij-cosijk*vkj[j]/lenkj)*4.184070*100.0;

      f_a[ii][j] += f1;
      f_a[kk][j] += f2;
      f_a[jj][j] += -f1-f2;
    }
  }

  return p_a_t;
}

double GOLMAA_PROTEINS2008_debug_ff_calcDIHE(double *crd,int numatom,double *p_d,double **f_d,double Kd1,double Kd2,double *dih_equ,int **pairs,int numdihe){
  int i,j,k,l;
  int ii,jj,kk,ll;

  double fi[3],fj[3],fk[3],fl[3];
  double m[3],n[3],m_n[3],n_n[3],lm,ln;
  double vij[3],vkj[3],vkl[3];
  double lkj;
  double vijvkj,vklvkj;

  double dvdpsi;
  
  double atom[4][3];
  double dihed;
  double p_d_t=0.0;

  double pi;

  pi=acos(-1.0);

  for (i=0;i<numdihe;++i) p_d[i] = 0.0;
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) f_d[i][j] = 0.0;
  
  for (i=0;i<numdihe;++i) {
    ii=pairs[i][0];
    jj=pairs[i][1];
    kk=pairs[i][2];
    ll=pairs[i][3];

    if ( ii < jj && jj < kk && kk < ll  ) {
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

      dihed=acos(inprod(m_n,n_n,3));
      if (inprod(vij,n,3)>0) dihed=-dihed;
      if (dihed<0.0) dihed=2.0*pi+dihed;
      
      vijvkj=inprod(vij,vkj,3);
      vklvkj=inprod(vkl,vkj,3);

      dvdpsi=Kd1*sin(dihed-dih_equ[i])+3.0*Kd2*sin(3.0*(dihed-dih_equ[i]));
      
      p_d[i] = Kd1*(1.0-cos(dihed-dih_equ[i]))+Kd2*(1.0-cos(3.0*(dihed-dih_equ[i])));
      p_d_t += p_d[i];

      for (j=0;j<3;++j) {
	fi[j] = -dvdpsi*lkj*m[j]/(lm*lm);
	fl[j] =  dvdpsi*lkj*n[j]/(ln*ln);
	fj[j] = (-fi[j]+(vijvkj/(lkj*lkj))*fi[j]-(vklvkj/(lkj*lkj))*fl[j]);
	fk[j] = (-fl[j]-(vijvkj/(lkj*lkj))*fi[j]+(vklvkj/(lkj*lkj))*fl[j]);
	
	f_d[ii][j] += fi[j]*4.184070*100.0;
	f_d[jj][j] += fj[j]*4.184070*100.0;
	f_d[kk][j] += fk[j]*4.184070*100.0;
	f_d[ll][j] += fl[j]*4.184070*100.0;
      }
    }
  }

  return p_d_t;

}

double GOLMAA_PROTEINS2008_debug_ff_calcff_nonlocal(double *crd, int numatom,int *NC_index,int numNC,int *NotNC_index,int numNotNC,
					      double *ALJ_natatt,double *BLJ_natatt,double e_natatt,double ALJ_Repul,
					      double *p_natatt,double *p_repul,double **f_natatt,double **f_repul) {
  int i,j,k;
  int atomi,atomj;

  double vec[3];
  double len,len2,len6,len12;
  double p12,p6;

  for (i=0;i<numNC;++i) p_natatt[i] = 0.0;

  for (i=0;i<numNotNC;++i) p_repul[i] = 0.0;

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      f_natatt[i][j] = 0.0;
      f_repul[i][j] = 0.0;
    }
  }

  for (i=0;i<numNC;++i) {
    atomi=NC_index[i*2+0];
    atomj=NC_index[i*2+1];

    len2 = 0.0;
    for(k=0;k<3;++k){
      vec[k] = crd[atomi*3+k]-crd[atomj*3+k];
      len2 += vec[k]*vec[k];
    }
    len = sqrt(len2);
    len6=len2;
    len12=len2;
    for (k=0;k<2;++k)  len6 = len6*len2;
    for (k=0;k<5;++k)  len12 = len12*len2;
    p12 = ALJ_natatt[i]/len12;
    p6 = BLJ_natatt[i]/len6;
    //    p_natatt[atomi] += e_natatt*(p12-2.0*p6);
    p_natatt[i] = e_natatt*(p12-2.0*p6);
    for (k=0;k<3;++k) {
      f_natatt[atomi][k] +=  e_natatt*(12.0*p12-6.0*2.0*p6)/(len2)*vec[k]*UNIT;
      f_natatt[atomj][k] += -e_natatt*(12.0*p12-6.0*2.0*p6)/(len2)*vec[k]*UNIT;
    }
  }

  for (i=0;i<numNotNC;++i) {
    atomi=NotNC_index[i*2+0];
    atomj=NotNC_index[i*2+1];

    len2 = 0.0;
    for(k=0;k<3;++k){
      vec[k] = crd[atomi*3+k]-crd[atomj*3+k];
      len2 += vec[k]*vec[k];
    }
    len = sqrt(len2);
    len12=len2;
    for (k=0;k<5;++k) len12 = len12*len2;
    p12 = ALJ_Repul/len12;
    //    p_repul[atomi] += p12;
    p_repul[i] = p12;
    for (k=0;k<3;++k) {
      f_repul[atomi][k] +=  12.0*p12/(len2)*vec[k]*UNIT;
      f_repul[atomj][k] += -12.0*p12/(len2)*vec[k]*UNIT;
    }
  }
}

