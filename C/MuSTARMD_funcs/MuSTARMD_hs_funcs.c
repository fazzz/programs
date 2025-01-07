
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "FFLc.h"
#include "FFL.h"
#include "MB.h"

#include "PTLb.h"

#include "LA.h"
#include "TOPO.h"
#include "mymath.h"
#include "EF.h"

#include "MuSTARMD_hs_funcs.h"

double ffLc_calcffandforce_HS(double *crd, int numatom,struct potential *ene,struct force *f,struct AmberParmL ap) {
  int i;
  int numnb,num14;
  double *n_d;

  numnb=(*ene).parm.numnb;
  num14=(*ene).parm.num14;

  (*ene).p_e_t=ffLc_calcffandforce_nb_HS((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,
					 (*f).f_e,numnb,
					 (*ene).parm.indexnb,numatom,crd);

  (*ene).p_e_14_t=ffLc_calcffandforce_14_HS((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,
					    (*f).f_e_14,num14,
					    (*ene).parm.index14,numatom,crd);

  (*ene).p_d_t=ffLc_calcDIHE_ffandforce_Cartesian_HS(ap.PH,ap.PA,
						     ap.PN,ap.PK,ap.PHASE,
						     (*f).f_d,crd,numatom,ap.NPHIH,ap.MPHIA);

  (*ene).p_a_t=ffLc_calcANGLE_ffandforce_Cartesian_HS(ap.TH,ap.TA,
						      ap.TK,ap.TEQ,
						      (*f).f_a,crd,
						      numatom,ap.NTHETH,ap.MTHETA);

  (*ene).p_b_t=ffLc_calcBOND_ffandforce_Cartesian_HS(ap.BH,ap.BA,
						     ap.RK,ap.REQ,
						     (*f).f_b,crd,numatom,ap.NBONH,ap.MBONA);

  (*ene).p_t=(*ene).p_e_t+(*ene).p_e_14_t+(*ene).p_d_t+(*ene).p_a_t+(*ene).p_b_t;

  return (*ene).p_t;
}

double ffLc_calcffandforce_nb_HS(double *ele, double *ALJ, double *BLJ, double *f,
				 int numnb, int *indexnb,
				 int num_atom,double *cord){
  int i;
  int num_a_prot,NUM_A_PROT;
  double len,len2,len6;
  double vec[3],fdummy[3];
  double pes,pt;
  double p12,p6;

  int i_3,i_3_1,i_3_2,i_2;
  int num_a_prot_3,num_a_prot_3_1,num_a_prot_3_2;
  int NUM_A_PROT_3,NUM_A_PROT_3_1,NUM_A_PROT_3_2;
  int num_a_prot_num_a_prot_NUM_A_PROT;

  double fes,fLJ,fes_LJ;

  double p_t=0.0;

  for(i=0;i<num_atom;++i) {
    i_3=i*3;
    i_3_1=i*3+1;
    i_3_2=i*3+2;

    f[i_3]=0.0;
    f[i_3_1]=0.0;
    f[i_3_2]=0.0;
  }

  for(i=0;i<numnb;++i){
    i_2=i*2;
    num_a_prot=indexnb[i_2];
    NUM_A_PROT=indexnb[i_2+1];

    num_a_prot_3=num_a_prot*3;
    NUM_A_PROT_3=NUM_A_PROT*3;
    num_a_prot_3_1=num_a_prot_3+1;
    NUM_A_PROT_3_1=NUM_A_PROT_3+1;
    num_a_prot_3_2=num_a_prot_3+2;
    NUM_A_PROT_3_2=NUM_A_PROT_3+2;

    vec[0] = cord[NUM_A_PROT_3]-cord[num_a_prot_3];
    vec[1] = cord[NUM_A_PROT_3_1]-cord[num_a_prot_3_1];
    vec[2] = cord[NUM_A_PROT_3_2]-cord[num_a_prot_3_2];
    len2 = vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
    len = sqrt(len2);
    len6 = len2*len2*len2;

    pes=ele[num_a_prot]*ele[NUM_A_PROT]/(len);

    num_a_prot_num_a_prot_NUM_A_PROT=num_a_prot*num_a_prot+NUM_A_PROT;
    p12 = ALJ[num_a_prot_num_a_prot_NUM_A_PROT]/(len6*len6);
    p6  = BLJ[num_a_prot_num_a_prot_NUM_A_PROT]/len6;

    pt=pes+p12-p6;

    p_t += pt;

    fes=-pes;
    fLJ=-6.0*(2.0*p12-p6);
    fes_LJ=(fLJ-fes)/len2*UNIT;

    fdummy[0] = fes_LJ*vec[0];
    fdummy[1] = fes_LJ*vec[1];
    fdummy[2] = fes_LJ*vec[2];

    f[num_a_prot_3]   +=fdummy[0];
    f[num_a_prot_3_1] +=fdummy[1];
    f[num_a_prot_3_2] +=fdummy[2];

    f[NUM_A_PROT_3]   -=fdummy[0];
    f[NUM_A_PROT_3_1] -=fdummy[1];
    f[NUM_A_PROT_3_2] -=fdummy[2];
  }

  return p_t;
}

double ffLc_calcffandforce_14_HS(double *ele, double *ALJ, double *BLJ, double *f,
				 int numnb, int *indexnb,
				 int num_atom,double *cord){
  int i;
  int num_a_prot,NUM_A_PROT;
  double len,len2,len6;
  double vec[3],fdummy[3];
  double pes,pt;
  double p12,p6;

  int i_3,i_3_1,i_3_2,i_2;
  int num_a_prot_3,num_a_prot_3_1,num_a_prot_3_2;
  int NUM_A_PROT_3,NUM_A_PROT_3_1,NUM_A_PROT_3_2;
  int num_a_prot_num_a_prot_NUM_A_PROT;

  double fes,fLJ,fes_LJ;

  double p_t=0.0;

 for(i=0;i<num_atom;++i) {
    i_3=i*3;
    i_3_1=i*3+1;
    i_3_2=i*3+2;

    f[i_3]=0.0;
    f[i_3_1]=0.0;
    f[i_3_2]=0.0;
  }

  for(i=0;i<numnb;++i){
    i_2=i*2;
    num_a_prot=indexnb[i_2];
    NUM_A_PROT=indexnb[i_2+1];

    num_a_prot_3=num_a_prot*3;
    NUM_A_PROT_3=NUM_A_PROT*3;
    num_a_prot_3_1=num_a_prot_3+1;
    NUM_A_PROT_3_1=NUM_A_PROT_3+1;
    num_a_prot_3_2=num_a_prot_3+2;
    NUM_A_PROT_3_2=NUM_A_PROT_3+2;

    vec[0] = cord[NUM_A_PROT_3]-cord[num_a_prot_3];
    vec[1] = cord[NUM_A_PROT_3_1]-cord[num_a_prot_3_1];
    vec[2] = cord[NUM_A_PROT_3_2]-cord[num_a_prot_3_2];
    len2 = vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
    len = sqrt(len2);
    len6 = len2*len2*len2;

    pes=1.0/1.2*ele[num_a_prot]*ele[NUM_A_PROT]/(len);

    num_a_prot_num_a_prot_NUM_A_PROT=num_a_prot*num_a_prot+NUM_A_PROT;
    p12 = 0.5*ALJ[num_a_prot_num_a_prot_NUM_A_PROT]/(len6*len6);
    p6  = 0.5*BLJ[num_a_prot_num_a_prot_NUM_A_PROT]/len6;

    pt=pes+p12-p6;

    p_t += pt;

    fes=-pes;
    fLJ=-6.0*(2.0*p12-p6);
    fes_LJ=(fLJ-fes)/len2*UNIT;

    fdummy[0] = fes_LJ*vec[0];
    fdummy[1] = fes_LJ*vec[1];
    fdummy[2] = fes_LJ*vec[2];

    f[num_a_prot_3]   +=fdummy[0];
    f[num_a_prot_3_1] +=fdummy[1];
    f[num_a_prot_3_2] +=fdummy[2];

    f[NUM_A_PROT_3]   -=fdummy[0];
    f[NUM_A_PROT_3_1] -=fdummy[1];
    f[NUM_A_PROT_3_2] -=fdummy[2];
  }

  return p_t;
}

double ffLc_calcDIHE_ffandforce_Cartesian_HS(int **PH, int **PA,
					     double *PN, double *PK, double *PHASE,
					     double *f_d,double *cord,
					     int numatom, int NPHIH, int MPHIA) {
  int i,j,k,l,flag;
  int dtype;

  double fa,fb[3],fc[3];
  double *n1,*n2,ln1,ln2;
  double vij[3],vkj[3],vki[3],vjl[3],vji[3],vik[3],vkl[3],vlk[3];
  double op1[3],op2[3],op3[3],op4[3],op5[3],op6[3];
  
  double atom[4][3];
  double cosdih,sindih;
  double dihedang;

  double cosat,sinat;

  double p_t;

  int na1,na2,na3,na4;
  int N;
  double K,phase;

  n1=(double *)gcemalloc(sizeof(double)*3);
  n2=(double *)gcemalloc(sizeof(double)*3);

  for (i=0;i<numatom*3;++i) f_d[i] = 0.0;
  
  for (i=0;i<NPHIH;++i) {
    dtype = PH[i][4]-1;
    na1=abs(PH[i][0]);
    na2=abs(PH[i][1]);
    na3=abs(PH[i][2]);
    na4=abs(PH[i][3]);

    K=PK[dtype];
    N=PN[dtype];
    phase=PHASE[dtype];

    for (k=0;k<3;++k) {
      atom[0][k]=cord[na1+k];
      atom[1][k]=cord[na2+k];
      atom[2][k]=cord[na3+k];
      atom[3][k]=cord[na4+k];
    }
    
    for (j=0;j<3;++j) {
      vij[j] = atom[1][j]-atom[0][j];
      vkj[j] = atom[1][j]-atom[2][j];
      vki[j] = atom[0][j]-atom[2][j];
      vjl[j] = atom[3][j]-atom[1][j];
      vji[j] = atom[0][j]-atom[1][j];
      vik[j] = atom[2][j]-atom[0][j];
      vkl[j] = atom[3][j]-atom[2][j];
      vlk[j] = atom[2][j]-atom[3][j];
    }

    outprod(vij,vkj,n1);
    outprod(vkj,vkl,n2);
    ln1=sqrt(inprod(n1,n1,3));
    ln2=sqrt(inprod(n2,n2,3));
  
    csdih(atom[0],atom[1],atom[2],atom[3],&cosdih,&sindih);

    if (N==1) {
      fa=-K*N*cos(phase)*UNIT;
      cosat=cosdih; sinat=sindih;
    }
    else if (N==2) {
      fa=-K*N*2.0*cosdih*cos(phase)*UNIT;
      cosat=cosdih*cosdih-sindih*sindih;
      sinat=2.0*sindih*cosdih;
    }
    else if (N==3) {
      fa=-K*N*(-4.0*sindih*sindih+3.0)*cos(phase)*UNIT;
      cosat=cosdih*cosdih*cosdih-3.0*sindih*sindih*cosdih;
      sinat=4.0*sindih*cosdih*cosdih;
    }
    else if (PN[dtype]==4) {
      fa=-K*N*4.0*(cosdih*(2.0*cosdih*cosdih-1.0))*cos(phase)*UNIT;
      cosat=cosdih*cosdih*cosdih*cosdih-7.0*sindih*sindih*cosdih*cosdih;
      sinat=8.0*sindih*cosdih*cosdih*cosdih;
    }
    else {
      printf("error:phase must be 1~4\n");
      exit(1);
    }

    p_t += K*(1.0+(cosat*cos(phase)+sinat*sin(phase)));

    for (j=0;j<3;++j) fb[j]=(n2[j]/ln2-cosdih*n1[j]/ln1)/ln1;
    for (j=0;j<3;++j) fc[j]=(n1[j]/ln1-cosdih*n2[j]/ln2)/ln2;

    outprod(fb,vkj,op1);
    outprod(fb,vki,op2);
    outprod(fc,vlk,op3);
    outprod(fb,vij,op4);
    outprod(fc,vjl,op5);
    outprod(fc,vkj,op6);

    for (j=0;j<3;++j) {
      f_d[na1+j] += fa*op1[j];
      f_d[na2+j] += fa*(-op2[j]+op3[j]);
      f_d[na3+j] += fa*(-op4[j]+op5[j]);
      f_d[na4+j] += fa*op6[j];
    }
  }
  
  for (i=0;i<MPHIA;++i) {
    dtype = PA[i][4]-1;
    na1=abs(PA[i][0]);
    na2=abs(PA[i][1]);
    na3=abs(PA[i][2]);
    na4=abs(PA[i][3]);

    K=PK[dtype];
    N=PN[dtype];
    phase=PHASE[dtype];

    for (k=0;k<3;++k) {
      atom[0][k]=cord[na1+k];
      atom[1][k]=cord[na2+k];
      atom[2][k]=cord[na3+k];
      atom[3][k]=cord[na4+k];
    }

    for (j=0;j<3;++j) {
      vij[j] = atom[1][j]-atom[0][j];
      vkj[j] = atom[1][j]-atom[2][j];
      vki[j] = atom[0][j]-atom[2][j];
      vjl[j] = atom[3][j]-atom[1][j];
      vji[j] = atom[0][j]-atom[1][j];
      vik[j] = atom[2][j]-atom[0][j];
      vkl[j] = atom[3][j]-atom[2][j];
      vlk[j] = atom[2][j]-atom[3][j];
    }

    outprod(vij,vkj,n1);
    outprod(vkj,vkl,n2);
    ln1=sqrt(inprod(n1,n1,3));
    ln2=sqrt(inprod(n2,n2,3));
  
    csdih(atom[0],atom[1],atom[2],atom[3],&cosdih,&sindih);

    if (N==1) {
      fa=-K*N*cos(phase)*UNIT;
      cosat=cosdih; sinat=sindih;
    }
    else if (N==2) {
      fa=-K*N*2.0*cosdih*cos(phase)*UNIT;
      cosat=cosdih*cosdih-sindih*sindih;
      sinat=2.0*sindih*cosdih;
    }
    else if (N==3) {
      fa=-K*N*(-4.0*sindih*sindih+3.0)*cos(phase);
      cosat=cosdih*cosdih*cosdih-3.0*sindih*sindih*cosdih;
      sinat=4.0*sindih*cosdih*cosdih;
    }
    else if (N==4) {
      fa=-K*N*4.0*(cosdih*(2.0*cosdih*cosdih-1.0))*cos(phase);
      cosat=cosdih*cosdih*cosdih*cosdih-7.0*sindih*sindih*cosdih*cosdih;
      sinat=8.0*sindih*cosdih*cosdih*cosdih;
    }
    else {
      printf("error:periodicity must be 1~4\n");
      exit(1);
    }

    p_t += K*(1.0+(cosat*cos(phase)+sinat*sin(phase)));

    for (j=0;j<3;++j) fb[j]=(n2[j]/ln2-cosdih*n1[j]/ln1)/ln1;
    for (j=0;j<3;++j) fc[j]=(n1[j]/ln1-cosdih*n2[j]/ln2)/ln2;

    outprod(fb,vkj,op1);
    outprod(fb,vki,op2);
    outprod(fc,vlk,op3);
    outprod(fb,vij,op4);
    outprod(fc,vjl,op5);
    outprod(fc,vkj,op6);

    for (j=0;j<3;++j) {
      f_d[na1+j] += fa*op1[j];
      f_d[na2+j] += fa*(-op2[j]+op3[j]);
      f_d[na3+j] += fa*(-op4[j]+op5[j]);
      f_d[na4+j] += fa*op6[j];
    }
  }

  return p_t;
}

double ffLc_calcANGLE_ffandforce_Cartesian_HS(int **TH, int **TA, 
					      double *TK, double *TEQ,
					      double *f_a,double *cord,
					      int numatom, int NTHETH, int MTHETA){
  int i,j,k,l;
  int type;
  double kang,ang_eq;
  double atom[3][3];

  double lenij,lenkj;
  double vij[3],vkj[3];
  double cosijk,angijk;
  double f1,f2;

  double p_t=0.0;

  int na1,na2,na3;

  for (i=0;i<numatom*3;++i) f_a[i] = 0.0;

  for (i=0;i<NTHETH;++i) {
    na1=abs(TH[i][0]);
    na2=abs(TH[i][1]);
    na3=abs(TH[i][2]);

    type = TH[i][3]-1;
    kang = TK[type];
    ang_eq = TEQ[type];

    for (k=0;k<3;++k) {
      atom[0][k]=cord[na1+k];
      atom[1][k]=cord[na2+k];
      atom[2][k]=cord[na3+k];
    }

    lenij=0.0;
    lenkj=0.0;
    for (j=0;j<3;++j) {
      vij[j]=atom[1][j]-atom[0][j];
      vkj[j]=atom[1][j]-atom[2][j];
      lenij += vij[j]*vij[j];
      lenkj += vkj[j]*vkj[j];
    }
    lenij=sqrt(lenij);
    lenkj=sqrt(lenkj);

    cosijk=inprod(vij,vkj,3);
    cosijk=cosijk/lenij/lenkj;
    angijk = acos(cosijk);

    p_t += kang*(angijk-ang_eq)*(angijk-ang_eq);

    for (j=0;j<3;++j) {
      f1 = -2.0*kang*(angijk-ang_eq)/(lenij*sin(angijk))*(vkj[j]/lenkj-cosijk*vij[j]/lenij)*UNIT;
      f2 = -2.0*kang*(angijk-ang_eq)/(lenkj*sin(angijk))*(vij[j]/lenij-cosijk*vkj[j]/lenkj)*UNIT;
      
      f_a[na1+j] += f1;
      f_a[na3+j] += f2;
      f_a[na2+j] += -f1-f2;
    }
  }

  for (i=0;i<MTHETA;++i) {
    na1=abs(TA[i][0]);
    na2=abs(TA[i][1]);
    na3=abs(TA[i][2]);

    type = TA[i][3]-1;
    kang = TK[type];
    ang_eq = TEQ[type];

    for (k=0;k<3;++k) {
      atom[0][k]=cord[na1+k];
      atom[1][k]=cord[na2+k];
      atom[2][k]=cord[na3+k];
    }

    lenij=0.0;
    lenkj=0.0;
    for (j=0;j<3;++j) {
      vij[j]=atom[1][j]-atom[0][j];
      vkj[j]=atom[1][j]-atom[2][j];
      lenij += vij[j]*vij[j];
      lenkj += vkj[j]*vkj[j];
    }
    lenij=sqrt(lenij);
    lenkj=sqrt(lenkj);

    cosijk=inprod(vij,vkj,3);
    cosijk=cosijk/lenij/lenkj;
    angijk = acos(cosijk);

    p_t += kang*(angijk-ang_eq)*(angijk-ang_eq);

    for (j=0;j<3;++j) {
      f1 = -2.0*kang*(angijk-ang_eq)/(lenij*sin(angijk))*(vkj[j]/lenkj-cosijk*vij[j]/lenij)*UNIT;
      f2 = -2.0*kang*(angijk-ang_eq)/(lenkj*sin(angijk))*(vij[j]/lenij-cosijk*vkj[j]/lenkj)*UNIT;
      
      f_a[na1+j] += f1;
      f_a[na3+j] += f2;
      f_a[na2+j] += -f1-f2;
    }
  }

  return 0;
}

double ffLc_calcBOND_ffandforce_Cartesian_HS(int **BH, int **BA,double *RK, double *REQ,
					     double *f_b,double *cord,
					     int numatom, int NBONH, int MBONA){
  int i,j,k;
  int type;
  double f;
  double lenij;
  double vecij[3];
  double atom[2][3];

  int na1,na2;

  double K,EQ;
  double p_t=0.0;

  for (i=0;i<numatom*3;++i) f_b[i] = 0.0;
  
  for (i=0;i<NBONH;++i) {
    type = BH[i][2]-1;
    na1=abs(BH[i][0]);
    na2=abs(BH[i][1]);

    K=RK[type];
    EQ=REQ[type];

    for (k=0;k<3;++k) {
      atom[0][k]=cord[na1+k];
      atom[1][k]=cord[na2+k];
    }
  
    lenij=0.0;
    for (j=0;j<3;++j) {
      vecij[j]=atom[1][j]-atom[0][j];
      lenij+=vecij[j]*vecij[j];
    }
    lenij=sqrt(lenij);

    p_t += K*(lenij-EQ)*(lenij-EQ);

    for (j=0;j<3;++j) {
      f = -2.0*K*(lenij-EQ)*vecij[j]/lenij*UNIT;
      f_b[na1+j] += f;
      f_b[na2+j] += -f;
    }
  }

  for (i=0;i<MBONA;++i) {
    type = BA[i][2]-1;
    na1=abs(BH[i][0]);
    na2=abs(BH[i][1]);

    K=RK[type];
    EQ=REQ[type];

    for (k=0;k<3;++k){
      atom[j][k]=cord[na1+k];
      atom[j][k]=cord[na2+k];
    }
  
    lenij=0.0;
    for (j=0;j<3;++j) {
      vecij[j]=atom[1][j]-atom[0][j];
      lenij+=vecij[j]*vecij[j];
    }
    lenij=sqrt(lenij);

    p_t += K*(lenij-EQ)*(lenij-EQ);

    for (j=0;j<3;++j) {
      f = -2.0*K*(lenij-EQ)*vecij[j]/lenij*UNIT;
      f_b[na1+j] += f;
      f_b[na2+j] += -f;
    }
  }

  return 0;
}

int ffLc_set_calcffandforce_HS(struct potential *ene, struct force *f,struct AmberParmL ap){
  int i,numatom,fd;
  int numnb,num14;
  int num_a_prot,NUM_A_PROT;
  char dummy;

  ffLc_set_non_bonding_index_1(&numnb,&num14,ap);
  (*ene).parm.numnb=numnb;
  (*ene).parm.num14=num14;
  (*ene).parm.indexnb=(int *)gcemalloc(sizeof(int)*numnb*2);
  (*ene).parm.index14=(int *)gcemalloc(sizeof(int)*num14*2);
  ffLc_set_non_bonding_index_2((*ene).parm.indexnb,(*ene).parm.index14,ap);

  numatom=ap.NATOM;
  (*ene).parm.ele=(double *)gcemalloc(sizeof(double)*numatom);
  (*ene).parm.ALJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);
  (*ene).parm.BLJ=(double *)gcemalloc(sizeof(double)*numatom*numatom);

  ffLc_set_NB_PARM((*ene).parm.ele,(*ene).parm.ALJ,(*ene).parm.BLJ,numatom,ap);

  (*f).f_t=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_e=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_e_14=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_d=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_a=(double *)gcemalloc(sizeof(double)*numatom*3);
  (*f).f_b=(double *)gcemalloc(sizeof(double)*numatom*3);

}
