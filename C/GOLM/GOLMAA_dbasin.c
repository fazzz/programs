
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLMAA.h"
#include "GOLMAA_set.h"
#include "GOLMAA_dbasin.h"
#include "PTL.h"

#define ON 1
#define OFF 0

double GOLMAAff_dbasin_calcff(double *crd, int numatom,struct potential_GOLMAA_dbasin *ene,double delta, double deltaV) {
  int i,j;
  double p_t;
  double p_1,p_2;
  double p_sub1,p_sub2,p_sub3,p_sub4;

  GOLMAAff_calcff(crd,numatom,&(ene->L1),OFF,ON,ON,ene->nb_matrix1);
  GOLMAAff_calcff(crd,numatom,&(ene->L2),OFF,ON,ON,ene->nb_matrix2);

  p_1=(ene->L1).p_t;
  p_2=(ene->L2).p_t;

  p_sub1=0.5*(p_1+p_2+deltaV);
  p_sub2=0.5*(p_1-p_2-deltaV);
  p_sub3=sqrt(p_sub2*p_sub2+delta*delta);
  p_sub4=p_1-p_2-deltaV;

  ene->p_t=p_sub1-p_sub3;

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      ene->f_t[i][j]=0.5*((ene->L1).f_t[i][j]+(ene->L2).f_t[i][j])
	-0.25*(p_sub4)/p_sub3*((ene->L1).f_t[i][j]-(ene->L2).f_t[i][j]);
    }
  }
  
  return p_t;
}

double GOLMAAff_dbasin_calcff_ratio(double *crd, int numatom,struct potential_GOLMAA_dbasin *ene,double delta, double deltaV) {
  int i,j;
  double p_1,p_2;
  double p_sub1,p_sub2,p_sub3,p_sub4;
  double ratio;

  GOLMAAff_calcff(crd,numatom,&(ene->L1),OFF,ON,ON,ene->nb_matrix1);
  GOLMAAff_calcff(crd,numatom,&(ene->L2),OFF,ON,ON,ene->nb_matrix2);

  p_1=(ene->L1).p_t;
  p_2=(ene->L2).p_t;

  p_sub1=0.5*(p_1+p_2+deltaV);
  p_sub2=0.5*(p_1-p_2-deltaV);
  p_sub3=sqrt(p_sub2*p_sub2+delta*delta);
  p_sub4=p_1-p_2-deltaV;

  ene->p_t=p_sub1-p_sub3;

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      ene->f_t[i][j]=0.5*((ene->L1).f_t[i][j]+(ene->L2).f_t[i][j])
	-0.25*(p_sub4)/p_sub3*((ene->L1).f_t[i][j]-(ene->L2).f_t[i][j]);
    }
  }

  ratio=log(-delta/(-p_1+ene->p_t));
  
  return ratio;
}


double GOLMAA_dbasin_calcTorque(double *Q,double *crd,struct potential_GOLMAA_dbasin *ene,int numclut,int *nNumClutOfParent,int *terminal,int *origin,double delta, double deltaV) {
  int i;
  double *Q1,*Q2;
  double pot_d_L1,pot_d_L2;
  double p_t;
  double p_1,p_2;
  double p_sub1,p_sub2,p_sub3,p_sub4;

  Q1=(double *)gcemalloc(sizeof(double)*numclut);
  Q2=(double *)gcemalloc(sizeof(double)*numclut);
  pot_d_L1=GOLMAA_calcTorque(Q1,crd,(ene->L1).DEQ,(ene->L1).FC_dihed,numclut,nNumClutOfParent,terminal,origin);
  pot_d_L2=GOLMAA_calcTorque(Q2,crd,(ene->L2).DEQ,(ene->L2).FC_dihed,numclut,nNumClutOfParent,terminal,origin);

  p_1=pot_d_L1;
  p_2=pot_d_L2;

  p_sub1=0.5*(p_1+p_2+deltaV);
  p_sub2=0.5*(p_1-p_2-deltaV);
  p_sub3=sqrt(p_sub2*p_sub2+delta*delta);
  p_sub4=p_1-p_2-deltaV;

  p_t=p_sub1-p_sub3;

  for (i=0;i<numclut;++i) Q[i]=0.5*(Q1[i]+Q2[i])-0.25*(p_sub4)/p_sub3*(Q1[i]-Q2[i]);

  return p_t;
}

int GOLMAAff_dbasin_set_calcff(struct potential_GOLMAA_dbasin *ene, double *refcrd1, double *refcrd2,int numatom,double R_C_D,double constant){
  int i;

  (*ene).nb_matrix1=(int **)gcemalloc(sizeof(int *)*numatom);
  (*ene).nb_matrix2=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) {
    (*ene).nb_matrix1[i]=(int *)gcemalloc(sizeof(int)*numatom);
    (*ene).nb_matrix2[i]=(int *)gcemalloc(sizeof(int)*numatom);
  }
  GOLMAAff_set_calcff(&((*ene).L1),refcrd1,numatom,(*ene).nb_matrix1,R_C_D,constant);
  GOLMAAff_set_calcff(&((*ene).L2),refcrd2,numatom,(*ene).nb_matrix2,R_C_D,constant);

  (*ene).f_t=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    (*ene).f_t[i]=(double *)gcemalloc(sizeof(double)*3);
  }

}

