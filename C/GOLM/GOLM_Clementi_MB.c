
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "GOLM_Clementi_set.h"
#include "GOLM_Clementi.h"
#include "GOLM_Clementi_MB.h"

double GOLM_Clementi_MB_ff_calcff(double *crd, int numatom, double de, double d2,
				  struct potential_GOLM_Clementi_MB *ene) {
  int i,j;
  double e1,e2;
  double **f1,**f2;
  double A,B,C,D;

  GOLM_Clementi_ff_calcff(crd,numatom,&((*ene).e1));
  GOLM_Clementi_ff_calcff(crd,numatom,&((*ene).e2));

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

double GOLM_Clementi_MB_ff_set_calcff(struct potential_GOLM_Clementi_MB *ene,
				      double *refcrd1,double *refcrd2,
				      double *refcrdAA1,double *refcrdAA2,
				      int numCAatom, int numatom, double ep) {
  int i,j;

  GOLM_Clementi_ff_set_calcff2(&((*ene).e1),refcrd1,refcrdAA1,numCAatom,numatom,ep);
  GOLM_Clementi_ff_set_calcff2(&((*ene).e2),refcrd2,refcrdAA2,numCAatom,numatom,ep);

  (*ene).f_MB = (double **)gcemalloc(sizeof(double *)*numCAatom);
  for (i=0;i<numCAatom;++i) (*ene).f_MB[i]=(double *)gcemalloc(sizeof(double)*3);

  return 0.0;
}

double GOLM_Clementi_MB_Kai(double *crd,int numatom, double de,double d, double d2, 
			    struct potential_GOLM_Clementi_MB *ene) {
  double e1,e2;
  double A,B,C,D;
  double Kai;

  GOLM_Clementi_ff_calcff(crd,numatom,&((*ene).e1));
  GOLM_Clementi_ff_calcff(crd,numatom,&((*ene).e2));

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

double MD_Propagetor_NH_MP1998_GOLM_Clementi_MB(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLM_Clementi_MB *e_GOLM, double de, double d2) {
  int i,j,k;
  double KE;

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  //  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*(*e_GOLM).e1.f_t[i][j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*(*e_GOLM).f_MB[i][j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  //  GOLM_Clementi_ff_calcff(crd,numatom,&((*e_GOLM).e1));
  GOLM_Clementi_MB_ff_calcff(crd,numatom,de,d2,e_GOLM);
  //  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*(*e_GOLM).e1.f_t[i][j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*(*e_GOLM).f_MB[i][j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}
