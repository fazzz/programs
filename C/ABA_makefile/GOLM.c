
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLM.h"
#include "NC.h"
#include "PTL.h"
#include "TOPO.h"
#include "MB.h"

double GOLMff_calcff(double *crd, int numatom,struct potential_GOLM *ene,
		     int flagb,int flaga, int flagd, int flagnc, int flagnn) {
  int i,j;

  if (flagb==ON) {
    GOLMpote_calcBOND((*ene).p_b,crd,(*ene).FC_bond,(*ene).BEQ,(*ene).num_bond,(*ene).index_bond);
    GOLMforc_calcBOND((*ene).f_b,crd,(*ene).FC_bond,(*ene).BEQ,(*ene).num_bond,(*ene).index_bond,numatom);
  }

  if (flaga==ON) {
    GOLMpote_calcANGLE((*ene).p_a,crd,(*ene).FC_angle,(*ene).AEQ,(*ene).num_angl,(*ene).index_angl);
    GOLMforc_calcANGLE((*ene).f_a,crd,(*ene).FC_angle,(*ene).AEQ,(*ene).num_angl,(*ene).index_angl,numatom);
  }

  if (flagd==ON) {
    GOLMpote_calcDIHE((*ene).p_d,crd,(*ene).DEQ,(*ene).FC_dihed1,(*ene).FC_dihed2,(*ene).num_dihe,(*ene).index_dihe);
    GOLMforc_calcDIHE((*ene).f_d,crd,(*ene).DEQ,(*ene).FC_dihed1,(*ene).FC_dihed2,(*ene).num_dihe,(*ene).index_dihe,numatom);
  }

  if (flagnc==ON) {
    GOLMpote_calcNatAtt((*ene).p_natatt,crd,(*ene).ALJ_natatt,(*ene).BLJ_natatt,(*ene).ep_natatt,(*ene).num_natatt,numatom,(*ene).index_natatt);
    GOLMforc_calcNatAtt((*ene).f_natatt,crd,(*ene).ALJ_natatt,(*ene).BLJ_natatt,(*ene).ep_natatt,(*ene).num_natatt,numatom,(*ene).index_natatt);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////

  if (flagnn==ON) {
    GOLMpote_calcRepul((*ene).p_repul,crd,(*ene).ALJ_repul,numatom);
    GOLMforc_calcRepul((*ene).f_natatt,crd,(*ene).ALJ_repul,numatom);
  }

  (*ene).p_t=0.0;
  (*ene).p_b_t=0.0;
  (*ene).p_a_t=0.0;
  (*ene).p_d_t=0.0;
  (*ene).p_natatt_t=0.0;
  (*ene).p_repul_t=0.0;

  for (i=0;i<(*ene).num_dihe;++i) {
    (*ene).p_t+=(*ene).p_d[i];
    (*ene).p_d_t+=(*ene).p_d[i];
  }
  for (i=0;i<(*ene).num_angl;++i) {
    (*ene).p_t+=(*ene).p_a[i];
    (*ene).p_a_t+=(*ene).p_a[i];
  }
  for (i=0;i<(*ene).num_bond;++i) {
    (*ene).p_t+=(*ene).p_b[i];
    (*ene).p_b_t+=(*ene).p_b[i];
  }
  for (i=0;i<numatom;++i) {
    (*ene).p_t+=(*ene).p_natatt[i];
    (*ene).p_natatt_t+=(*ene).p_natatt[i];
  }
  for (i=0;i<numatom;++i) {
    (*ene).p_t+=(*ene).p_repul[i];
    (*ene).p_repul_t+=(*ene).p_repul[i];
  }

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      (*ene).f_t[i][j]=(*ene).f_b[i][j]+(*ene).f_a[i][j]+(*ene).f_d[i][j]+(*ene).f_natatt[i][j]+(*ene).f_repul[i][j];
    }
  }
  
  return (*ene).p_t;
}

int GOLMpote_calcBOND(double *p_b,double *cord, double FC_bond, double *BEQ,int num_bond,int *index_bond){
  int i,j,k;
  double length;
  double atom[2][3];

  for (i=0;i<num_bond;++i) {
    for (j=0;j<2;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[(index_bond[i*2+j])*3+k];

    length = len(atom[0],atom[1]);
    p_b[i] = FC_bond*(length-BEQ[i])*(length-BEQ[i]);
  }

  return 0;
}

void GOLMforc_calcBOND(double **f_b,double *cord, double FC_bond, double *BEQ,int num_bond,int *index_bond, int numatom){
  int i,j,k;
  double f;
  double length;
  double atom[2][3];
  
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) f_b[i][j] = 0.0;
  for (i=0;i<num_bond;++i) {
    for (j=0;j<2;++j)
      for (k=0;k<3;++k)
	atom[j][k]=cord[index_bond[i*2+j]*3+k];

    length = len(atom[0],atom[1]);
    for (j=0;j<3;++j) {
      f = -2.0*FC_bond*(length-BEQ[i])*(atom[1][j]-atom[0][j])/length*UNIT;
      f_b[(index_bond[i*2])][j] += -f;
      f_b[(index_bond[i*2+1])][j] += f;
    }
  }
}

int GOLMpote_calcANGLE(double *p_a,double *crd,double FC_angle,double *AEQ,int num_angl,int *index_angl){
  int i,j,k;
  double atom[3][3];
  double angle,dang;
  double pi;
  
  pi=acos(-1.0);
  
  for (i=0;i<num_angl;++i) p_a[i] = 0.0;
  
  for (i=0;i<num_angl;++i) {
     for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	atom[j][k]=crd[(index_angl[i*3+j])*3+k];
  
     angle = ang(atom[0],atom[1],atom[2]);
     if (angle>pi) angle-=2.0*pi;
     else if (angle<-1.0*pi) angle+=2.0*pi;

     if ((dang=angle-AEQ[i])>pi) dang-=2.0*pi;
     else if ((dang=angle-AEQ[i])<-1.0*pi) dang+=2.0*pi;

     p_a[i] = FC_angle*dang*dang;
  }

  return 0;
}

int GOLMforc_calcANGLE(double **f_a,double *crd,double FC_angle,double *AEQ,int num_angl,int *index_angl,int numatom){
  int i,j,k;
  double atom[3][3];
  double vij[3],vkj[3];
  double lenij,lenkj;
  double cosijk,angijk;
  double f1,f2;
  
  for (i=0;i<numatom;++i) 
    for (j=0;j<3;++j)
      f_a[i][j] = 0.0;  

  for (i=0;i<num_angl;++i) {
    for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	atom[j][k]=crd[(index_angl[i*3+j])*3+k];

    lenij = len(atom[0],atom[1]);
    lenkj = len(atom[2],atom[1]);
    for (j=0;j<3;++j) {
      vij[j]=atom[1][j]-atom[0][j];
      vkj[j]=atom[1][j]-atom[2][j];
    }
    cosijk=inprod(vij,vkj,3);
    cosijk=cosijk/lenij/lenkj;
    angijk = acos(cosijk);

    for (j=0;j<3;++j) {
      f1 = -2.0*FC_angle*(angijk-AEQ[i])/(lenij*sin(angijk))*(vkj[j]/lenkj-cosijk*vij[j]/lenij)*UNIT;
      f2 = -2.0*FC_angle*(angijk-AEQ[i])/(lenkj*sin(angijk))*(vij[j]/lenij-cosijk*vkj[j]/lenkj)*UNIT;

      f_a[(index_angl[i*3])][j] += f1;
      f_a[(index_angl[i*3+2])][j] += f2;
      f_a[(index_angl[i*3+1])][j] += -f1-f2;
    }
  }

  return 0;
}

void GOLMpote_calcDIHE(double *p_d,double *crd,double *DEQ,double FC_dihed1,double FC_dihed2,int num_dihed,int *index_dihed){
  int i,j,k,l;

  double atom[4][3];
  double cosdih,sindih;
  double dihedang,delta;
  double pi;

  pi=acos(-1.0);
  for (i=0;i<num_dihed;++i)
    p_d[i]=0.0;

  for (i=0;i<num_dihed;++i) {
    for (j=0;j<4;++j)
      for (k=0;k<3;++k)
	atom[j][k]=crd[index_dihed[i*4+j]*3+k];
  
    dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);
    if (dihedang>pi) dihedang-=2.0*pi;
    else if (dihedang<-1.0*pi) dihedang+=2.0*pi;

    if ((delta=dihedang-DEQ[i])>pi) delta-=2.0*pi;
    else if ((delta=dihedang-DEQ[i])<-1.0*pi) delta+=2.0*pi;

    p_d[i] = FC_dihed1*(1.0-cos(delta));
  }
}

void GOLMforc_calcDIHE(double **f_d,double *crd,double *DEQ,double FC_dihed1,double FC_dihed2,int num_dihed,int *index_dihed,int numatom) {
  int i,j,k,l;

  double fa,fb[3],fc[3];
  double *n1,*n2,ln1,ln2;
  double vij[3],vkj[3],vki[3],vjl[3],vji[3],vik[3],vkl[3];
  double op1[3],op2[3],op3[3],op4[3],op5[3],op6[3];
  
  double atom[4][3];
  double cosdih,sindih;
  double dihedang;

  n1=(double *)gcemalloc(sizeof(double)*3);
  n2=(double *)gcemalloc(sizeof(double)*3);

  for (i=0;i<numatom;++i) 
    for (j=0;j<3;++j)
      f_d[i][j] = 0.0;
  for (i=0;i<num_dihed;++i) {
    for (j=0;j<4;++j)
      for (k=0;k<3;++k)
	atom[j][k]=crd[index_dihed[i*4+j]*3+k];
  
    dihedang = pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0.0);

    for (j=0;j<3;++j) {
      vij[j] = atom[1][j]-atom[0][j];
      vkj[j] = atom[1][j]-atom[2][j];
      vki[j] = atom[0][j]-atom[2][j];
      vjl[j] = atom[3][j]-atom[1][j];
      vji[j] = atom[0][j]-atom[1][j];
      vik[j] = atom[2][j]-atom[0][j];
      vkl[j] = atom[3][j]-atom[2][j];
    }

    outprod(vij,vkj,n1);
    outprod(vkj,vkl,n2);
    ln1=sqrt(inprod(n1,n1,3));
    ln2=sqrt(inprod(n2,n2,3));
  
    csdih(atom[0],atom[1],atom[2],atom[3],&cosdih,&sindih);

    fa=-FC_dihed2*cos(DEQ[i]);

    for (j=0;j<3;++j) fb[j]=(n2[j]/ln2-cosdih*n1[j]/ln1)/ln1;
    for (j=0;j<3;++j) fc[j]=(n1[j]/ln1-cosdih*n2[j]/ln2)/ln2;

    outprod(fb,vkj,op1);
    outprod(fc,vki,op2);
    outprod(fb,vik,op3);
    outprod(fb,vij,op4);
    outprod(fc,vjl,op5);
    outprod(fc,vkj,op6);

    for (j=0;j<3;++j) {
      f_d[(index_dihed[i*4])][j] += fa*op1[j]*UNIT;
      f_d[(index_dihed[i*4+1])][j] += fa*(-op2[j]+op3[j])*UNIT;
      f_d[(index_dihed[i*4+2])][j] += fa*(-op4[j]+op5[j])*UNIT;
      f_d[(index_dihed[i*4+3])][j] += fa*op6[j]*UNIT;
    }
  }
}

int GOLMpote_calcNatAtt(double *p_natatt,double *crd,double *ALJ_natatt,double *BLJ_natatt,double ep,int num_natatt,int numatom,int *index_natatt) {
  int i,j;
  int num_a_prot,NUM_A_PROT;
  double vec[3];
  double len,len2,len10,len12;
  double p12,p10;

  for(i=0;i<numatom;++i) p_natatt[i]=0.0;

  for(i=0;i<num_natatt;++i){
    num_a_prot=index_natatt[i*2];
    NUM_A_PROT=index_natatt[i*2+1];
    len2 = 0.0;
    for(j=0;j<3;++j){
      vec[j] = crd[NUM_A_PROT*3+j]-crd[num_a_prot*3+j];
      len2 += vec[j]*vec[j];
    }
    len = sqrt(len2);
    len10=len2;
    len12=len2;
    for (j=0;j<4;++j)  len10 = len10*len2;
    for (j=0;j<5;++j)  len12 = len12*len2;
    p12 = 5.0*ALJ_natatt[i]/len12;
    p10 = 6.0*BLJ_natatt[i]/len10;
    if (1.0<p12-p10+1.0) {
      p_natatt[num_a_prot] += ep*1.0;
      p_natatt[NUM_A_PROT] += ep*1.0;
    }
    else {
      p_natatt[num_a_prot] += ep*(p12-p10+1.0);
      p_natatt[NUM_A_PROT] += ep*(p12-p10+1.0);
    }
  }

  return 0;
}

int GOLMpote_calcRepul(double *p_repul,double *crd,double ALJ_repul,int numatom) {
  int i,j,k;
  int ii,jj;
  double vec[3];
  double len,len2,len12;
  double p12;

  for(i=0;i<numatom;++i) p_repul[i]=0.0;

  ii=0;
  for (i=0;i<numatom;++i) {
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      ++ii;
    }
    jj=ii;
    for (j=i;j<numatom;++j) {
      if ( strncmp(AP.IGRAPH[j],"CA",2)==0 ) {
  	++jj;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0 && strncmp(AP.IGRAPH[j],"CA",2)==0) {
  	if (jj>ii+3) {
	  len2 = 0.0;
	  for(k=0;k<3;++k){
	    vec[k] = crd[j*3+k]-crd[i*3+k];
	    len2 += vec[k]*vec[k];
	  }
	  len = sqrt(len2);
	  len12=len2;
	  for (k=0;k<5;++k)  len12 = len12*len2;
	  p12 = ALJ_repul/len12;
	  p_repul[i] += p12;
	  p_repul[j] += p12;
  	}
      }
    }
  }

  return 0;
}

int GOLMforc_calcNatAtt(double **f_natatt,double *crd,double *ALJ_natatt,double *BLJ_natatt,double ep,int num_natatt,int numatom,int *index_natatt) {
  int i,j;
  int num_a_prot,NUM_A_PROT;
  double vec[3];
  double len,len2,len10,len12;
  double p12,p10;

  for(i=0;i<numatom;++i) for(j=0;j<3;++j) f_natatt[i][j]=0.0;

  for(i=0;i<num_natatt;++i){
    num_a_prot=index_natatt[i*2];
    NUM_A_PROT=index_natatt[i*2+1];
    len2 = 0.0;
    for(j=0;j<3;++j){
      vec[j] = crd[NUM_A_PROT*3+j]-crd[num_a_prot*3+j];
      len2 += vec[j]*vec[j];
    }
    len = sqrt(len2);
    len10=len2;
    len12=len2;
    for (j=0;j<4;++j)  len10 = len10*len2;
    for (j=0;j<5;++j)  len12 = len12*len2;
    p12 = 5.0*ALJ_natatt[i]/len12;
    p10 = 6.0*BLJ_natatt[i]/len10;
    if (1.0<p12-p10+1.0) {
      for (j=0;j<3;++j) {
	f_natatt[num_a_prot][j] += -(12.0*p12-10.*p10)/(len2)*vec[j]*UNIT*ep;
	f_natatt[NUM_A_PROT][j] += (12.0*p12-10.*p10)/(len2)*vec[j]*UNIT*ep;
      }
    }
  }

  return 0;
}

int GOLMforc_calcRepul(double **f_repul,double *crd,double ALJ_repul,int numatom) {
  int i,j,k;
  int ii,jj;
  double vec[3];
  double len,len2,len12;
  double p12;

  for(i=0;i<numatom;++i) for(j=0;j<3;++j) f_repul[i][j]=0.0;

  ii=0;
  for (i=0;i<numatom;++i) {
    if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
      ++ii;
    }
    jj=ii;
    for (j=i;j<numatom;++j) {
      if ( strncmp(AP.IGRAPH[j],"CA",2)==0 ) {
  	++jj;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0 && strncmp(AP.IGRAPH[j],"CA",2)==0) {
  	if (jj>ii+3) {
	  len2 = 0.0;
	  for(k=0;k<3;++k){
	    vec[k] = crd[j*3+k]-crd[i*3+k];
	    len2 += vec[k]*vec[k];
	  }
	  len = sqrt(len2);
	  len12=len2;
	  for (k=0;k<5;++k)  len12 = len12*len2;
	  p12 = ALJ_repul/len12;
	  for (k=0;k<3;++k) {
	    f_repul[i][k] += -(12.0*p12)/(len2)*vec[k]*UNIT;
	    f_repul[j][k] += (12.0*p12)/(len2)*vec[k]*UNIT;
	  }
  	}
      }
    }
  }

  return 0;
}


