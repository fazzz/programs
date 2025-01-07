
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MD_NHC_MP1996.h"
#include "EF.h"

#define nys 3

double MD_Propagetor_NH_Single_set_MP1996(int nc,double dt,double *dt2,double wdt2[3],double wdt4[3]) {
  int i,j;

  double omega1=1.0/(2.0-1.0/pow(2.0,3));
  double omega2=1.0/(2.0-1.0/pow(2.0,3));
  double omega3=1.0-2.0/(2.0-1.0/pow(2.0,3));

  *dt2=0.5*dt;

  wdt2[0]=omega1*dt/2.0/nc;
  wdt2[1]=omega2*dt/2.0/nc;
  wdt2[2]=omega3*dt/2.0/nc;

  wdt4[0]=omega1*dt/4.0/nc;
  wdt4[1]=omega2*dt/4.0/nc;
  wdt4[2]=omega3*dt/4.0/nc;

}

double MD_Propagetor_NH_Single_part_MP1996(double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,int nc,double wdt4[3],double wdt2[3]) {
  int i,j;
  double Scale=1.0,AA,G,KE;

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  G=(2.0*KE-NfKT)/Q;

  for (i=0;i<nc;++i) {
    for (j=0;j<nys;++j) {
      *V_zeta+=G*wdt4[j];
      AA=exp(-wdt2[j]*(*V_zeta));
      Scale=Scale*AA;
      G=(Scale*Scale*2.0*KE-NfKT)/Q;
      *zeta+=(*V_zeta)*wdt2[j];
      *V_zeta+=G*wdt4[j];
    }
  }

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=Scale*vel[i*3+j];

}

double MD_Propagetor_NH_MP1998_GOLM_Clementi(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLM_Clementi *e_GOLM) {
  int i,j,k;
  double KE;

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*(*e_GOLM).f_t[i][j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  GOLM_Clementi_ff_calcff(crd,numatom,e_GOLM);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*(*e_GOLM).f_t[i][j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}

double MD_Propagetor_NH_MP1998_GOLMAA_JCTC2011(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f, struct potential_GOLMAA_hybrid *e_GOLM) {
  int i,j,k;
  double *frc,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  ffL_calcffandforce_14DAB_woH(crd,numatom,e,f);
  GOLMAA_hyb_ff_calcff(crd,numatom,e_GOLM);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}

/***********************************************************************************************************************************************************************************************************************************************************************/
/* double MD_Propagetor_NH_MP1998_14LJdab(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f) { */
/*   int i,j,k;																															       */
/*   double *frc,KE;																														       */
/* 																																       */
/*   frc=(double *)gcemalloc(sizeof(double)*numatom*3);																										       */
/* 																																       */
/*   for (i=0;i<numatom;++i)																													       */
/*     for (j=0;j<3;++j) 																													       */
/*       frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j];																					       */
/* 																																       */
/*   MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);																					       */
/* 																																       */
/*   //////////////////////////////////////////////////////////////////////////////////////																					       */
/*   for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];																						       */
/*   for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];																							       */
/*   ffL_calcffandforce_14DAB_woH(crd,numatom,e,f);																										       */
/*   for (i=0;i<numatom;++i)																													       */
/*     for (j=0;j<3;++j) 																													       */
/*       frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j];																					       */
/*   for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];																						       */
/*   //////////////////////////////////////////////////////////////////////////////////////																					       */
/* 																																       */
/*   MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);																					       */
/* 																																       */
/*   KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];																					       */
/* 																																       */
/*   *KEv=0.5*Q*(*V_zeta)*(*V_zeta);																												       */
/*   *PEv=NfKT*(*zeta);																														       */
/* 																																       */
/*   return KE;																															       */
/* }																																       */
/***********************************************************************************************************************************************************************************************************************************************************************/


double MD_Propagetor_NH_MP1998_GOLMAA_JCTC2011_wonat(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f, struct potential_GOLMAA_hybrid *e_GOLM) {
  int i,j,k;
  double *frc,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j]+(*e_GOLM).f_repul[i][j];

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  ffL_calcffandforce_14DAB_woH(crd,numatom,e,f);
  GOLMAA_hyb_ff_calcff(crd,numatom,e_GOLM);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j]+(*e_GOLM).f_repul[i][j];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}

double MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLMAA_PROTEINS2008 *e_GOLM) {
  int i,j,k;
  double *frc,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_GOLM).f_b[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_d[i][j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  GOLMAA_PROTEINS2008_ff_calcff(crd,numatom,e_GOLM);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_GOLM).f_b[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_d[i][j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}

/*******************************************************************************************************************************************************************************************************************************************************************************************/
/* double MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLMAA_PROTEINS2008 *e_GOLM) { */
/*   int i,j,k;																																		   */
/*   double *frc,KE;																																	   */
/* 																																			   */
/*   frc=(double *)gcemalloc(sizeof(double)*numatom*3);																													   */
/* 																																			   */
/*   for (i=0;i<numatom;++i)																																   */
/*     for (j=0;j<3;++j) 																																   */
/*       frc[i*3+j]=(*e_GOLM).f_b[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_d[i][j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];																			   */
/* 																																			   */
/*   MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);																								   */
/* 																																			   */
/*   //////////////////////////////////////////////////////////////////////////////////////																								   */
/*   for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];																									   */
/*   for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];																										   */
/*   GOLMAA_PROTEINS2008_ff_calcff(crd,numatom,e_GOLM);																													   */
/*   for (i=0;i<numatom;++i)																																   */
/*     for (j=0;j<3;++j) 																																   */
/*       frc[i*3+j]=(*e_GOLM).f_b[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_d[i][j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];																			   */
/*   for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];																									   */
/*   //////////////////////////////////////////////////////////////////////////////////////																								   */
/* 																																			   */
/*   MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);																								   */
/* 																																			   */
/*   KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];																								   */
/* 																																			   */
/*   *KEv=0.5*Q*(*V_zeta)*(*V_zeta);																															   */
/*   *PEv=NfKT*(*zeta);																																	   */
/* 																																			   */
/*   return KE;																																		   */
/* }																																			   */
/*******************************************************************************************************************************************************************************************************************************************************************************************/


double MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008_b(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLMAA_PROTEINS2008 *e_GOLM) {
  int i,j,k;
  double *frc,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_GOLM).f_b[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_d[i][j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  GOLMAA_PROTEINS2008_ff_calcff_b(crd,numatom,e_GOLM);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_GOLM).f_b[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_d[i][j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}

double MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008_b_2013_08_31(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLMAA_PROTEINS2008 *e_GOLM,double fact_b,double fact_a,double fact_t,double fact_NC,double fact_NNC) {
  int i,j,k;
  double *frc,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=fact_b*(*e_GOLM).f_b[i][j]+fact_a*(*e_GOLM).f_a[i][j]+fact_t*(*e_GOLM).f_d[i][j]+fact_NC*(*e_GOLM).f_natatt[i][j]+fact_NNC*(*e_GOLM).f_repul[i][j];

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  GOLMAA_PROTEINS2008_ff_calcff_b(crd,numatom,e_GOLM);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=fact_b*(*e_GOLM).f_b[i][j]+fact_a*(*e_GOLM).f_a[i][j]+fact_t*(*e_GOLM).f_d[i][j]+fact_NC*(*e_GOLM).f_natatt[i][j]+fact_NNC*(*e_GOLM).f_repul[i][j];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}

double MD_Propagetor_NH_MP1998_GOLMAA_MB_PROTEINS2008(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3],double de, double d2, struct potential_GOLMAA_MB_PROTEINS2008 *e_GOLM) {
  int i,j,k;
  double *frc,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) frc[i*3+j]=(*e_GOLM).f_MB[i][j];

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  GOLMAA_MB_PROTEINS2008_ff_calcff(crd,numatom,de,d2,e_GOLM);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) frc[i*3+j]=(*e_GOLM).f_MB[i][j];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}


double MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008_debug(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLMAA_PROTEINS2008 *e_GOLM,int bflag, int aflag, int dflag, int nflag) {
  int i,j,k;
  double *frc,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_GOLM).f_b[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_d[i][j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  GOLMAA_PROTEINS2008_debug_ff_calcff(crd,numatom,e_GOLM,bflag,aflag,dflag,nflag,bflag,aflag,dflag,nflag);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_GOLM).f_b[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_d[i][j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}

double MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008_debug2(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLMAA_PROTEINS2008 *e_GOLM,int numstep) {
  int i,j,k;
  double *frc,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_GOLM).f_b[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_d[i][j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  GOLMAA_PROTEINS2008_ff_calcff_debug(crd,numatom,e_GOLM,numstep);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_GOLM).f_b[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_d[i][j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}


double MD_Propagetor_NH_MP1998_14LJdab(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f) {
  int i,j,k;
  double *frc,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j];

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  ffL_calcffandforce_14DAB_woH(crd,numatom,e,f);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}

double MD_Propagetor_NH_MP1998_AAFF_Amber(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f) {
  int i,j,k;
  double *frc,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_e[i*3+j]+(*f).f_LJ[i*3+j]+(*f).f_e_14[i*3+j]+(*f).f_LJ_14[i*3+j];

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  ffL_calcffandforce(crd,numatom,e,f);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_e[i*3+j]+(*f).f_LJ[i*3+j]+(*f).f_e_14[i*3+j]+(*f).f_LJ_14[i*3+j];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}

double MD_Propagetor_NH_MP1998_AAFF_Amber_wflag(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f,int bodflag,int angflag,int dihflag,int LJ14flag,int es14flag,int LJflag,int esflag) {
  int i,j,k;
  double *frc,KE;

  int ON=1;
  int OFF=0;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) frc[i*3+j]=0.0;

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      if (esflag==ON)   frc[i*3+j]+=(*f).f_e[i*3+j];
      if (LJflag==ON)   frc[i*3+j]+=(*f).f_LJ[i*3+j];
      if (es14flag==ON) frc[i*3+j]+=(*f).f_e_14[i*3+j];
      if (LJ14flag==ON) frc[i*3+j]+=(*f).f_LJ_14[i*3+j];
      if (dihflag==ON)  frc[i*3+j]+=(*f).f_d[i*3+j];
      if (angflag==ON)  frc[i*3+j]+=(*f).f_a[i*3+j];
      if (bodflag==ON)  frc[i*3+j]-=(*f).f_b[i*3+j];
    }
  }

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  ffL_calcffandforce(crd,numatom,e,f);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      if (esflag==ON)   frc[i*3+j]+=(*f).f_e[i*3+j];
      if (LJflag==ON)   frc[i*3+j]+=(*f).f_LJ[i*3+j];
      if (es14flag==ON) frc[i*3+j]+=(*f).f_e_14[i*3+j];
      if (LJ14flag==ON) frc[i*3+j]+=(*f).f_LJ_14[i*3+j];
      if (dihflag==ON)  frc[i*3+j]+=(*f).f_d[i*3+j];
      if (angflag==ON)  frc[i*3+j]+=(*f).f_a[i*3+j];
      if (bodflag==ON)  frc[i*3+j]-=(*f).f_b[i*3+j];
    }
  }
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}

double MD_Propagetor_NH_MP1998_AAFF_Amber_UMB_wflag(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f,int bodflag,int angflag,int dihflag,int LJ14flag,int es14flag,int LJflag,int esflag,int *pairp,int nump,double *fcp,double *dih_equ,double **fUMB) {
  int i,j,k;
  double *frc,KE;
  double *pUMB;

  int ON=1;
  int OFF=0;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) frc[i*3+j]=0.0;

  pUMB=(double *)gcemalloc(sizeof(double)*nump);

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      if (esflag==ON)   frc[i*3+j]+=(*f).f_e[i*3+j];
      if (LJflag==ON)   frc[i*3+j]+=(*f).f_LJ[i*3+j];
      if (es14flag==ON) frc[i*3+j]+=(*f).f_e_14[i*3+j];
      if (LJ14flag==ON) frc[i*3+j]+=(*f).f_LJ_14[i*3+j];
      if (dihflag==ON)  frc[i*3+j]+=(*f).f_d[i*3+j];
      if (angflag==ON)  frc[i*3+j]+=(*f).f_a[i*3+j];
      if (bodflag==ON)  frc[i*3+j]-=(*f).f_b[i*3+j];
      frc[i*3+j]+=fUMB[i][j];
    }
  }

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  ffL_calcffandforce(crd,numatom,e,f);
  UMB_calc_dihetype_ff(crd,numatom,pairp,nump,fcp,dih_equ,pUMB,fUMB); 
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      if (esflag==ON)   frc[i*3+j]+=(*f).f_e[i*3+j];
      if (LJflag==ON)   frc[i*3+j]+=(*f).f_LJ[i*3+j];
      if (es14flag==ON) frc[i*3+j]+=(*f).f_e_14[i*3+j];
      if (LJ14flag==ON) frc[i*3+j]+=(*f).f_LJ_14[i*3+j];
      if (dihflag==ON)  frc[i*3+j]+=(*f).f_d[i*3+j];
      if (angflag==ON)  frc[i*3+j]+=(*f).f_a[i*3+j];
      if (bodflag==ON)  frc[i*3+j]-=(*f).f_b[i*3+j];
      frc[i*3+j]+=fUMB[i][j];
    }
  }
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}

double MD_Propagetor_NH_MP1998_AAFF_Amber_UMB(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f,int *pairp,int nump,double *fcp,double *dih_equ,double *pUMB,double **fUMB) {
  int i,j,k;
  double *frc,KE;

  int ON=1;
  int OFF=0;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) frc[i*3+j]=0.0;

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j]+(*f).f_e_14[i*3+j]+(*f).f_LJ[i*3+j]+(*f).f_e[i*3+j]+fUMB[i][j];
    }
  }

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  ffL_calcffandforce(crd,numatom,e,f);
  UMB_calc_dihetype_ff(crd,numatom,pairp,nump,fcp,dih_equ,pUMB,fUMB); 
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j]+(*f).f_e_14[i*3+j]+(*f).f_LJ[i*3+j]+(*f).f_e[i*3+j]+fUMB[i][j];
    }
  }
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}


double MD_Propagetor_NH_MP1998_AAFF_Amber_woH_wflag(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f,int bodflag,int angflag,int dihflag,int LJ14flag,int es14flag,int LJflag,int esflag) {
  int i,j,k;
  double *frc,KE;

  int ON=1;
  int OFF=0;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) frc[i*3+j]=0.0;

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      if (esflag==ON)   frc[i*3+j]+=(*f).f_e[i*3+j];
      if (LJflag==ON)   frc[i*3+j]+=(*f).f_LJ[i*3+j];
      if (es14flag==ON) frc[i*3+j]+=(*f).f_e_14[i*3+j];
      if (LJ14flag==ON) frc[i*3+j]+=(*f).f_LJ_14[i*3+j];
      if (dihflag==ON)  frc[i*3+j]+=(*f).f_d[i*3+j];
      if (angflag==ON)  frc[i*3+j]+=(*f).f_a[i*3+j];
      if (bodflag==ON)  frc[i*3+j]-=(*f).f_b[i*3+j];
    }
  }

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  ffL_calcffandforce_woH(crd,numatom,e,f);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      if (esflag==ON)   frc[i*3+j]+=(*f).f_e[i*3+j];
      if (LJflag==ON)   frc[i*3+j]+=(*f).f_LJ[i*3+j];
      if (es14flag==ON) frc[i*3+j]+=(*f).f_e_14[i*3+j];
      if (LJ14flag==ON) frc[i*3+j]+=(*f).f_LJ_14[i*3+j];
      if (dihflag==ON)  frc[i*3+j]+=(*f).f_d[i*3+j];
      if (angflag==ON)  frc[i*3+j]+=(*f).f_a[i*3+j];
      if (bodflag==ON)  frc[i*3+j]-=(*f).f_b[i*3+j];
    }
  }
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}

double MD_Propagetor_NH_MP1998_AAFF_Amber_woH_UMB_wflag(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f,int bodflag,int angflag,int dihflag,int LJ14flag,int es14flag,int LJflag,int esflag,int *pairp,int nump,double *fcp,double *dih_equ,double **fUMB) {
  int i,j,k;
  double *frc,KE;
  double *pUMB;

  int ON=1;
  int OFF=0;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) frc[i*3+j]=0.0;

  pUMB=(double *)gcemalloc(sizeof(double)*nump);

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      if (esflag==ON)   frc[i*3+j]+=(*f).f_e[i*3+j];
      if (LJflag==ON)   frc[i*3+j]+=(*f).f_LJ[i*3+j];
      if (es14flag==ON) frc[i*3+j]+=(*f).f_e_14[i*3+j];
      if (LJ14flag==ON) frc[i*3+j]+=(*f).f_LJ_14[i*3+j];
      if (dihflag==ON)  frc[i*3+j]+=(*f).f_d[i*3+j];
      if (angflag==ON)  frc[i*3+j]+=(*f).f_a[i*3+j];
      if (bodflag==ON)  frc[i*3+j]-=(*f).f_b[i*3+j];
      frc[i*3+j]+=fUMB[i][j];
    }
  }

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  ffL_calcffandforce_woH(crd,numatom,e,f);
  UMB_calc_dihetype_ff(crd,numatom,pairp,nump,fcp,dih_equ,pUMB,fUMB); 
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      if (esflag==ON)   frc[i*3+j]+=(*f).f_e[i*3+j];
      if (LJflag==ON)   frc[i*3+j]+=(*f).f_LJ[i*3+j];
      if (es14flag==ON) frc[i*3+j]+=(*f).f_e_14[i*3+j];
      if (LJ14flag==ON) frc[i*3+j]+=(*f).f_LJ_14[i*3+j];
      if (dihflag==ON)  frc[i*3+j]+=(*f).f_d[i*3+j];
      if (angflag==ON)  frc[i*3+j]+=(*f).f_a[i*3+j];
      if (bodflag==ON)  frc[i*3+j]-=(*f).f_b[i*3+j];
      frc[i*3+j]+=fUMB[i][j];
    }
  }
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}

double MD_Propagetor_NH_MP1998_AAFF_AMBER(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e,struct force *f) {
  int i,j,k;
  double *frc,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_e_14[i*3+j]+(*f).f_LJ_14[i*3+j]+(*f).f_e[i*3+j]+(*f).f_LJ[i*3+j];

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  ffL_calcffandforce(crd,numatom,e,f);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_e_14[i*3+j]+(*f).f_LJ_14[i*3+j]+(*f).f_e[i*3+j]+(*f).f_LJ[i*3+j];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}

// 2012-03-22
/*double MD_Propagetor_NH_MP1998_MB_FASYS(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3],double de, double d2, struct potential_GOLMAA_PROTEINS2008 *e_GOLM,int *dih_equ1,int *dih_equ2,double *e1,double *e2,double *p_MB,double **f_MB) {
  int i,j,k;
  double *frc,KE;
  double A,B,C,D;
  double p_d1_t,p_d2_t,*p_d1,*p_d2,**f_d1,**f_d2,**f1,**f2;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

  p_d1=(double *)gcemalloc(sizeof(double)*(*e_GOLM).num_dihe);
  p_d2=(double *)gcemalloc(sizeof(double)*(*e_GOLM).num_dihe);

  f1=(double **)gcemalloc(sizeof(double *)*numatom);
  f2=(double **)gcemalloc(sizeof(double *)*numatom);
  f_d1=(double **)gcemalloc(sizeof(double *)*numatom);
  f_d2=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) {
    f1[i]=(double *)gcemalloc(sizeof(double)*3);
    f2[i]=(double *)gcemalloc(sizeof(double)*3);
    f_d1[i]=(double *)gcemalloc(sizeof(double)*3);
    f_d2[i]=(double *)gcemalloc(sizeof(double)*3);
  }

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) frc[i*3+j]=f_MB[i][j];

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  //////////////////////////////////////////////////////////////////////////////////////

  (*e_GOLM).p_b_t=GOLMAA_PROTEINS2008_ff_calcBOND(crd,numatom,(*e_GOLM).p_b,(*e_GOLM).f_b,(*e_GOLM).Kb,(*e_GOLM).bon_equ,(*e_GOLM).pairs_bond,(*e_GOLM).num_bond);

  (*e_GOLM).p_a_t=GOLMAA_PROTEINS2008_ff_calcANGLE(crd,numatom,(*e_GOLM).p_a,(*e_GOLM).f_a,(*e_GOLM).Ka,(*e_GOLM).ang_equ,(*e_GOLM).pairs_angl,(*e_GOLM).num_angl);


  p_d1_t=GOLMAA_PROTEINS2008_ff_calcDIHE(crd,numatom,p_d1,f_d1,
					 (*e_GOLM).Kd1,0.0,0.0,dih_equ1,(*e_GOLM).pairs_dihe,
					 (*e_GOLM).num_dihe,(*e_GOLM).impindex);

  p_d2_t=GOLMAA_PROTEINS2008_ff_calcDIHE(crd,numatom,p_d2,f_d2,
					 (*e_GOLM).Kd1,0.0,0.0,dih_equ2,(*e_GOLM).pairs_dihe,
					 (*e_GOLM).num_dihe,(*e_GOLM).impindex);

  (*e1)=p_d1_t+(*e_GOLM).p_a_t+(*e_GOLM).p_b_t;
  (*e2)=p_d2_t+(*e_GOLM).p_a_t+(*e_GOLM).p_b_t;

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      f1[i][j]=f_d1[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_b[i][j];
      f2[i][j]=f_d2[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_b[i][j];
    }
  }

  A=0.5*((*e1)+(*e2)+de);
  B=0.5*((*e1)-(*e2)-de);
  C=sqrt(B*B+d2);
  D=(*e1)-(*e2)-de;
  
  *p_MB=A-C;
  //*p_MB=e1;

  for (i=0;i<numatom;++i) 
    for (j=0;j<3;++j) 
      f_MB[i][j]=0.5*(f1[i][j]+f2[i][j])-0.25*D/C*(f1[i][j]-f2[i][j]);
      //      f_MB[i][j]=f1[i][j];

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) frc[i*3+j]=f_MB[i][j];
  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
  }*/

// 2012-03-22
//double MD_Propagetor_NH_MP1998_FASYS(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3],double de, double d2, struct potential_GOLMAA_PROTEINS2008 *e_GOLM,int *dih_equ1,int *dih_equ2,double *p_MB,double **f_MB,double *KE) {
//  int i,j,k;
//  double *frc;
//  double A,B,C,D;
//  double p_d1_t,p_d2_t,*p_d1,*p_d2,**f_d1,**f_d2,**f1,**f2,e1,e2;

//  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

//  p_d1=(double *)gcemalloc(sizeof(double)*(*e_GOLM).num_dihe);
//  p_d2=(double *)gcemalloc(sizeof(double)*(*e_GOLM).num_dihe);

//  f1=(double **)gcemalloc(sizeof(double *)*numatom);
//  f2=(double **)gcemalloc(sizeof(double *)*numatom);
//  f_d1=(double **)gcemalloc(sizeof(double *)*numatom);
//  f_d2=(double **)gcemalloc(sizeof(double *)*numatom);
//  for (i=0;i<numatom;++i) {
//    f1[i]=(double *)gcemalloc(sizeof(double)*3);
//    f2[i]=(double *)gcemalloc(sizeof(double)*3);
//    f_d1[i]=(double *)gcemalloc(sizeof(double)*3);
//    f_d2[i]=(double *)gcemalloc(sizeof(double)*3);
//  }

//  for (i=0;i<numatom;++i) for (j=0;j<3;++j) frc[i*3+j]=f_MB[i][j];

//  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
//  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
//  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  //////////////////////////////////////////////////////////////////////////////////////

//  (*e_GOLM).p_b_t=GOLMAA_PROTEINS2008_ff_calcBOND(crd,numatom,(*e_GOLM).p_b,(*e_GOLM).f_b,(*e_GOLM).Kb,(*e_GOLM).bon_equ,(*e_GOLM).pairs_bond,(*e_GOLM).num_bond);

//  (*e_GOLM).p_a_t=GOLMAA_PROTEINS2008_ff_calcANGLE(crd,numatom,(*e_GOLM).p_a,(*e_GOLM).f_a,(*e_GOLM).Ka,(*e_GOLM).ang_equ,(*e_GOLM).pairs_angl,(*e_GOLM).num_angl);


//  p_d1_t=GOLMAA_PROTEINS2008_ff_calcDIHE(crd,numatom,p_d1,f_d1,
//					 (*e_GOLM).Kd1,0.0,0.0,dih_equ1,(*e_GOLM).pairs_dihe,
//					 (*e_GOLM).num_dihe,(*e_GOLM).impindex);

//  p_d2_t=GOLMAA_PROTEINS2008_ff_calcDIHE(crd,numatom,p_d2,f_d2,
//					 (*e_GOLM).Kd1,0.0,0.0,dih_equ2,(*e_GOLM).pairs_dihe,
//					 (*e_GOLM).num_dihe,(*e_GOLM).impindex);

//  e1=p_d1_t+(*e_GOLM).p_a_t+(*e_GOLM).p_b_t;
//  e2=p_d2_t+(*e_GOLM).p_a_t+(*e_GOLM).p_b_t;

//  for (i=0;i<numatom;++i) {
//    for (j=0;j<3;++j) {
//      f1[i][j]=f_d1[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_b[i][j];
//      f2[i][j]=f_d2[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_b[i][j];
//    }
//  }

//  A=0.5*(e1+e2+de);
//  B=0.5*(e1-e2-de);
//  C=sqrt(B*B+d2);
//  D=e1-e2-de;
  
//  //  *p_MB=A-C;
//  *p_MB=e1;

//  for (i=0;i<numatom;++i) 
//    for (j=0;j<3;++j) 
//      //      f_MB[i][j]=0.5*(f1[i][j]+f2[i][j])-0.25*D/C*(f1[i][j]-f2[i][j]);
//      f_MB[i][j]=f1[i][j];

//  for (i=0;i<numatom;++i) for (j=0;j<3;++j) frc[i*3+j]=f_MB[i][j];
//  //////////////////////////////////////////////////////////////////////////////////////
//  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
//  //////////////////////////////////////////////////////////////////////////////////////

//  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

//  *KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) *KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

//  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
//  *PEv=NfKT*(*zeta);

//  return 0.0/*KE*/;
//}

//2012-03-22
//double MD_Propagetor_NH_MP1998_FASYS_2(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3],double de, double d2, struct potential_GOLMAA_PROTEINS2008 *e_GOLM,int *dih_equ1,int *dih_equ2,double *p_MB,double **f_MB) {
//  int i,j,k;
//  double *frc,KE;
//  double A,B,C,D;
// double p_d1_t,p_d2_t,*p_d1,*p_d2,**f_d1,**f_d2,**f1,**f2,e1,e2;

// frc=(double *)gcemalloc(sizeof(double)*numatom*3);

// p_d1=(double *)gcemalloc(sizeof(double)*(*e_GOLM).num_dihe);
// p_d2=(double *)gcemalloc(sizeof(double)*(*e_GOLM).num_dihe);

// f1=(double **)gcemalloc(sizeof(double *)*numatom);
// f2=(double **)gcemalloc(sizeof(double *)*numatom);
// f_d1=(double **)gcemalloc(sizeof(double *)*numatom);
//f_d2=(double **)gcemalloc(sizeof(double *)*numatom);
//for (i=0;i<numatom;++i) {
//  f1[i]=(double *)gcemalloc(sizeof(double)*3);
//  f2[i]=(double *)gcemalloc(sizeof(double)*3);
//  f_d1[i]=(double *)gcemalloc(sizeof(double)*3);
//  f_d2[i]=(double *)gcemalloc(sizeof(double)*3);
//}

//for (i=0;i<numatom;++i) for (j=0;j<3;++j) frc[i*3+j]=f_MB[i][j];

//MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
//  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
// for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  //////////////////////////////////////////////////////////////////////////////////////

//(*e_GOLM).p_b_t=GOLMAA_PROTEINS2008_ff_calcBOND(crd,numatom,(*e_GOLM).p_b,(*e_GOLM).f_b,(*e_GOLM).Kb,(*e_GOLM).bon_equ,(*e_GOLM).pairs_bond,(*e_GOLM).num_bond);

//  (*e_GOLM).p_a_t=GOLMAA_PROTEINS2008_ff_calcANGLE(crd,numatom,(*e_GOLM).p_a,(*e_GOLM).f_a,(*e_GOLM).Ka,(*e_GOLM).ang_equ,(*e_GOLM).pairs_angl,(*e_GOLM).num_angl);


//  p_d1_t=GOLMAA_PROTEINS2008_ff_calcDIHE(crd,numatom,p_d1,f_d1,
//					 (*e_GOLM).Kd1,0.0,0.0,dih_equ1,(*e_GOLM).pairs_dihe,
//					 (*e_GOLM).num_dihe,(*e_GOLM).impindex);

//  p_d2_t=GOLMAA_PROTEINS2008_ff_calcDIHE(crd,numatom,p_d2,f_d2,
//					 (*e_GOLM).Kd1,0.0,0.0,dih_equ2,(*e_GOLM).pairs_dihe,
//					 (*e_GOLM).num_dihe,(*e_GOLM).impindex);

//  e1=p_d1_t+(*e_GOLM).p_a_t+(*e_GOLM).p_b_t;
//e2=p_d2_t+(*e_GOLM).p_a_t+(*e_GOLM).p_b_t;

//for (i=0;i<numatom;++i) {
//  for (j=0;j<3;++j) {
//    f1[i][j]=f_d1[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_b[i][j];
//    f2[i][j]=f_d2[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_b[i][j];
//  }
//}

//A=0.5*(e1+e2+de);
//B=0.5*(e1-e2-de);
//C=sqrt(B*B+d2);
//D=e1-e2-de;
  
////  *p_MB=A-C;
//*p_MB=e2;

//for (i=0;i<numatom;++i) 
//  for (j=0;j<3;++j) 
//    //      f_MB[i][j]=0.5*(f1[i][j]+f2[i][j])-0.25*D/C*(f1[i][j]-f2[i][j]);
//    f_MB[i][j]=f2[i][j];

//  for (i=0;i<numatom;++i) for (j=0;j<3;++j) frc[i*3+j]=f_MB[i][j];
  //////////////////////////////////////////////////////////////////////////////////////
  //for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  ////////////////////////////////////////////////////////////////////////////////////////

//  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

//  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

//  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
//*PEv=NfKT*(*zeta);

// return KE;
//}
