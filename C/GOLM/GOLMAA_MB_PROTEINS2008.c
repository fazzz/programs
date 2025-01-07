
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLMAA_PROTEINS2008.h"
#include "GOLMAA_MB_PROTEINS2008.h"

double GOLMAA_MB_PROTEINS2008_ff_calcff(double *crd, int numatom, double de, double d2,
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

  Kai=log(d/(e2+de-(*ene).p_MB));

  return Kai;
}

double GOLMAA_MB_PROTEINS2008_ff_calcff_wobaimp(double *crd, int numatom, double de, double d2,
						struct potential_GOLMAA_MB_PROTEINS2008 *ene) {
  int i,j;
  double e1,e2;
  double **f1,**f2;
  double A,B,C,D;

  GOLMAA_PROTEINS2008_ff_calcff_wobaimp(crd,numatom,&((*ene).e1));
  GOLMAA_PROTEINS2008_ff_calcff_wobaimp(crd,numatom,&((*ene).e2));

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

double GOLMAA_MB_PROTEINS2008_wobaimp_Kai(double *crd,int numatom, double de,double d, double d2, struct potential_GOLMAA_MB_PROTEINS2008 *ene) {
  double e1,e2;
  double A,B,C,D;
  double Kai;

  GOLMAA_PROTEINS2008_ff_calcff_wobaimp(crd,numatom,&((*ene).e1));
  GOLMAA_PROTEINS2008_ff_calcff_wobaimp(crd,numatom,&((*ene).e2));

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

double GOLMAA_MB_PROTEINS2008_wobaimp_Kai_debug(double *crd,int numatom, double de,double d, double d2, struct potential_GOLMAA_MB_PROTEINS2008 *ene) {
  double e1,e2;
  double A,B,C,D;
  double Kai;

  GOLMAA_PROTEINS2008_ff_calcff_wobaimp(crd,numatom,&((*ene).e1));
  GOLMAA_PROTEINS2008_ff_calcff_wobaimp(crd,numatom,&((*ene).e2));

  //  e1=(*ene).e1.p_natatt_t;
  //  e2=(*ene).e2.p_natatt_t;
  
  //  e1=(*ene).e1.p_d_t;
  //  e2=(*ene).e2.p_d_t;

  e1=(*ene).e1.p_repul_t+(*ene).e1.p_d_t;
  e2=(*ene).e2.p_repul_t+(*ene).e2.p_d_t;

  A=0.5*(e1+e2+de);
  B=0.5*(e1-e2-de);
  C=sqrt(B*B+d2);
  D=e1-e2-de;
  
  (*ene).p_MB=A-C;

  Kai=log((e1-(*ene).p_MB)/d);

  return Kai;
}


double GOLMAA_MB_PROTEINS2008_ff_calcff_harmo(double x, 
					      double de, double d2,
					      double k1, double x1,
					      double k2, double x2,
					      double *f) {
  int i,j;
  double e1,e2,e;
  double *f1,*f2;
  double A,B,C,D;

  e1=0.5*k1*(x-x1)*(x-x1);
  e2=0.5*k2*(x-x2)*(x-x2);

  f1=(double *)gcemalloc(sizeof(double)*3);
  f2=(double *)gcemalloc(sizeof(double)*3);

  for (i=0;i<3;++i) {
    f1[i]=f1[i];
    f2[i]=f2[i];
  }

  A=0.5*(e1+e2+de);
  B=0.5*(e1-e2-de);
  C=sqrt(B*B+d2);
  D=e1-e2-de;
  
  e=A-C;

  for (i=0;i<3;++i) 
    f[i]=0.5*(f1[i]+f2[i])-0.25*D/C*(f1[i]-f2[i]);

  return e;
}

double GOLMAA_MB_PROTEINS2008_harmo_Kai(double x, 
					double de, double d,double d2,
					double k1, double x1,
					double k2, double x2) {
  double e1,e2,e;
  double A,B,C,D;
  double Kai;

  e1=0.5*k1*(x-x1)*(x-x1);
  e2=0.5*k2*(x-x2)*(x-x2);

  A=0.5*(e1+e2+de);
  B=0.5*(e1-e2-de);
  C=sqrt(B*B+d2);
  D=e1-e2-de;
  
  e=A-C;

  Kai=log((e1-e)/d);

  return Kai;
}
