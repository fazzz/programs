
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLM_Clementi.h"

#include "TOPO.h"
#include "LA.h"
#include "mymath.h"

#include "PTL.h"

double GOLM_Clementi_ff_calcff(double *crd, int numatom,struct potential_GOLM_Clementi *ene) {
  int i,j;

  (*ene).p_b_t=GOLM_Clementi_ff_calcBOND(crd,numatom,(*ene).p_b,(*ene).f_b,(*ene).Kb,(*ene).bon_equ);

  (*ene).p_a_t=GOLM_Clementi_ff_calcANGLE(crd,numatom,(*ene).p_a,(*ene).f_a,(*ene).Ka,(*ene).ang_equ);

  (*ene).p_d_t=GOLM_Clementi_ff_calcDIHE(crd,numatom,(*ene).p_d,(*ene).f_d,(*ene).Kd1,(*ene).Kd2,(*ene).dih_equ);

  (*ene).p_natatt_t=GOLM_Clementi_ff_calcff_natatt(crd,numatom,
						   (*ene).index_natatt,(*ene).num_natatt,
						   (*ene).ALJ_natatt,(*ene).BLJ_natatt,(*ene).ep_natatt,
						   (*ene).p_natatt,(*ene).f_natatt);

  (*ene).p_repul_t=GOLM_Clementi_ff_calcff_non_natatt(crd,numatom,
						      (*ene).index_repul,(*ene).num_repul,
						      (*ene).ALJ_repul,
						      (*ene).p_repul,(*ene).f_repul);
  
  (*ene).p_t=(*ene).p_b_t+(*ene).p_a_t+(*ene).p_d_t+(*ene).p_natatt_t+(*ene).p_repul_t;
  for (i=0;i<numatom;++i) 
    for (j=0;j<3;++j) 
      (*ene).f_t[i][j]=(*ene).f_b[i][j]+(*ene).f_a[i][j]+(*ene).f_d[i][j]+(*ene).f_natatt[i][j]+(*ene).f_repul[i][j];

  return (*ene).p_t;
}

double GOLM_Clementi_ff_calcff_natatt(double *crd, int numatom,int *index_numatt,int numatt,
				      double *ALJ_natatt,double *BLJ_natatt,double ep_natatt,
				      double *p_natatt,double **f_natatt) {
  int i,j,k;
  int atomi,atomj;

  double vec[3];
  double len,len2,len10,len12;
  double p12,p10,f;
  double p_natatt_t=0.0;

  for (i=0;i<numatt;++i) p_natatt[i] = 0.0;
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) f_natatt[i][j] = 0.0;

  for (i=0;i<numatt;++i) {
    atomi=index_numatt[i*2+0];
    atomj=index_numatt[i*2+1];

    len2 = 0.0;
    for(j=0;j<3;++j){
      vec[j] = crd[atomi*3+j]-crd[atomj*3+j];
      len2 += vec[j]*vec[j];
    }
    len = sqrt(len2);
    len10=len2;
    len12=len2;
    for (j=0;j<4;++j)  len10 = len10*len2;
    for (j=0;j<5;++j)  len12 = len12*len2;
    p12 = ALJ_natatt[i]/len12;
    p10 = BLJ_natatt[i]/len10;
    p_natatt[i] = ep_natatt*(5.0*p12-6.0*p10);
    p_natatt_t += p_natatt[i];
    for (j=0;j<3;++j) {
      f = ep_natatt*(12.0*5.0*p12-10.0*6.0*p10)/(len2)*vec[j]*4.184070*100.0;
      f_natatt[atomi][j] +=  f;
      f_natatt[atomj][j] += -f;
    }
  }

  return p_natatt_t;
}

double GOLM_Clementi_ff_calcff_non_natatt(double *crd, int numatom,
					  int *index_notnumatt,int numnotatt,
					  double ALJ_repul,
					  double *p_repul,double **f_repul) {
  int i,j,k;
  int atomi,atomj;

  double vec[3];
  double len,len2,len12;
  double p12,f;
  double p_repul_t=0.0;

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) f_repul[i][j]=0.0;

  for (i=0;i<numnotatt;++i) {
    atomi=index_notnumatt[i*2+0];
    atomj=index_notnumatt[i*2+1];

    len2 = 0.0;
    for(j=0;j<3;++j){
      vec[j] = crd[atomi*3+j]-crd[atomj*3+j];
      len2 += vec[j]*vec[j];
    }
    len = sqrt(len2);
    len12=len2;
    for (j=0;j<5;++j) len12 = len12*len2;
    p12 = ALJ_repul/len12;
    p_repul[i] = p12;
    p_repul_t += p12;
    for (j=0;j<3;++j) {
      f = 12.0*p12/(len2)*vec[j]*4.184070*100.0;
      f_repul[atomi][j] +=  f;
      f_repul[atomj][j] += -f;
    }
  }

  return p_repul_t;
}

double GOLM_Clementi_ff_calcBOND(double *crd,int numatom,double *p_b,double **f_b,double Kb,double *bon_equ){
  int i,j,k;
  double f;
  double lenij;
  double atom[2][3];
  double p_b_t=0.0;

  for (i=0;i<numatom-1;++i) p_b[i] = 0.0;
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) f_b[i][j] = 0.0;
  
  for (i=0;i<numatom-1;++i) {
    for (j=0;j<3;++j) {
      atom[0][j]=crd[i*3+j];
      atom[1][j]=crd[(i+1)*3+j];
    }
  
    lenij = len(atom[0],atom[1]);
    p_b[i]=Kb*(lenij-bon_equ[i])*(lenij-bon_equ[i]);
    p_b_t += p_b[i];
    for (j=0;j<3;++j) {
      f = 2.0*Kb*(lenij-bon_equ[i])*(atom[1][j]-atom[0][j])/lenij*4.184070*100.0;
      f_b[i][j] += f;
      f_b[i+1][j] += -f;
    }
  }

  return p_b_t;
}

double GOLM_Clementi_ff_calcANGLE(double *crd,int numatom,double *p_a,double **f_a,double Ka,double *ang_equ){
  int i,j,k,l;
  double atom[3][3];
  double lenij,lenkj;
  double vij[3],vkj[3];
  double cosijk,angijk;
  double f1,f2;
  double p_a_t=0.0;
  double pi;

  pi=acos(-1.0);

  for (i=0;i<numatom-2;++i) p_a[i] = 0.0;
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) f_a[i][j] = 0.0;

  for (i=0;i<numatom-2;++i) {
    for (j=0;j<3;++j) {
      atom[0][j]=crd[i*3+j];
      atom[1][j]=crd[(i+1)*3+j];
      atom[2][j]=crd[(i+2)*3+j];
    }

    lenij = len(atom[0],atom[1]);
    lenkj = len(atom[2],atom[1]);
    for (j=0;j<3;++j) {
      vij[j]=atom[1][j]-atom[0][j];
      vkj[j]=atom[1][j]-atom[2][j];
    }
    cosijk=inprod(vij,vkj,3);
    cosijk=cosijk/lenij/lenkj;
    if (cosijk>=1.0) angijk=0.0;
    else if (cosijk<=0.0) angijk=pi;
    else  angijk = acos(cosijk);

    angijk = ang(atom[0],atom[1],atom[2]);

    p_a[i] = Ka*(angijk-ang_equ[i])*(angijk-ang_equ[i]);
    p_a_t += p_a[i];
    
    for (j=0;j<3;++j) {
      f1 = -2.0*Ka*(angijk-ang_equ[i])/(lenij*sin(angijk))*(vkj[j]/lenkj-cosijk*vij[j]/lenij)*4.184070*100.0;
      f2 = -2.0*Ka*(angijk-ang_equ[i])/(lenkj*sin(angijk))*(vij[j]/lenij-cosijk*vkj[j]/lenkj)*4.184070*100.0;

      f_a[i][j] += f1;
      f_a[i+2][j] += f2;
      f_a[i+1][j] += -f1-f2;
    }
  }

  return p_a_t;
}

double GOLM_Clementi_ff_calcDIHE(double *crd,int numatom,double *p_d,double **f_d,double Kd1,double Kd2,double *dih_equ){
  int i,j,k,l;

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

  for (i=0;i<numatom-3;++i) p_d[i] = 0.0;
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) f_d[i][j] = 0.0;
  
  for (i=0;i<numatom-3;++i) {
    for (j=0;j<3;++j) {
      atom[0][j]=crd[i*3+j];
      atom[1][j]=crd[(i+1)*3+j];
      atom[2][j]=crd[(i+2)*3+j];
      atom[3][j]=crd[(i+3)*3+j];
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
    if (dihed>=1.0) dihed=0.0;
    else if (dihed<=0.0) dihed=pi;
    else dihed=acos(dihed);

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

      f_d[i][j] += fi[j]*4.184070*100.0;
      f_d[i+1][j] += fj[j]*4.184070*100.0;
      f_d[i+2][j] += fk[j]*4.184070*100.0;
      f_d[i+3][j] += fl[j]*4.184070*100.0;
    }
  }

  return p_d_t;

}
