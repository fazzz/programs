
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLMAA_PROTEINS2008.h"
#include "GOLMAA_MB_PROTEINS2008.h"

double GOLMAA_MB_PROTEINS2008_ff_calcff_harmo(double *crd, int numatom, double de, double d2,
					      struct potential_GOLMAA_MB_PROTEINS2008 *ene) {
  int i,j;
  double e1,e2;
  double **f1,**f2;
  double A,B,C,D;

  GOLMAA_PROTEINS2008_ff_calcff(crd,numatom,&((*ene).e1));
  GOLMAA_PROTEINS2008_ff_calcff(crd,numatom,&((*ene).e2));

  f1=(double **)gcemalloc(sizeof(double *)*numatom);
  f2=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    f1[i]=(double *)gcemalloc(sizeof(double)*3);
    f2[i]=(double *)gcemalloc(sizeof(double)*3);
  }

  e1=(*ene).e1.p_t;
  e2=(*ene).e2.p_t;

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      f1[i][j]=(*ene).e1.f_t[i][j];
      f2[i][j]=(*ene).e2.f_t[i][j];
    }
  }

  A=0.5*(e1+e2+de);
  B=0.5*(e1-e2-de);
  C=sqrt(B*B+d2);
  D=e1-e2-de;
  
  (*ene).p_MB=A-C;

  for (i=0;i<numatom;++i) 
    for (j=0;j<3;++j) 
      (*ene).f_MB[i][j]=0.5*(f1[i][j]+f2[i][j])-0.25*D/C*(f1[i][j]-f2[i][j]);

  return 0.0;
}

double GOLMAA_MB_PROTEINS2008_ff_calcff_set(struct potential_GOLMAA_MB_PROTEINS2008 *ene, double *refcrd1,double *refcrd2, int numatom,int numres,int *non_bonding_index, int numnonbond, double ep, int nibnum,double criteria) {
  int i,j;

  GOLMAA_PROTEINS2008_ff_set_calcff(&((*ene).e1),refcrd1,numatom,numres,non_bonding_index,numnonbond,ep,nibnum,criteria);
  GOLMAA_PROTEINS2008_ff_set_calcff(&((*ene).e2),refcrd2,numatom,numres,non_bonding_index,numnonbond,ep,nibnum,criteria);

  (*ene).f_MB = (double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) (*ene).f_MB[i]=(double *)gcemalloc(sizeof(double)*3);

  return 0.0;
}

double GOLMAA_MB_PROTEINS2008_Kai(double *crd,int numatom, double de,double d, double d2, struct potential_GOLMAA_MB_PROTEINS2008 *ene) {
  double e1,e2;
  double A,B,C,D;
  double Kai;

  GOLMAA_PROTEINS2008_ff_calcff(crd,numatom,&((*ene).e1));
  GOLMAA_PROTEINS2008_ff_calcff(crd,numatom,&((*ene).e2));

  e1=(*ene).e1.p_t;
  e2=(*ene).e2.p_t;

  A=0.5*(e1+e2+de);
  B=0.5*(e1-e2-de);
  C=sqrt(B*B+d2);
  D=e1-e2-de;
  
  (*ene).p_MB=A-C;

  Kai=log((e1-(*ene).p_MB)/d);

  return Kai;
}

