
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MD_NHC_MP1996.h"
#include "EF.h"

#define nys 3

double MD_Generate_inivelo(double *vel,double *mass,int numatom,double KbT) {
  int i,j;
  double KE=0.0;
  double UNITT=418.4070;
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  for (i=0;i<numatom;++i) 
    for (j=0;j<3;++j) 
      vel[i*3+j]=Box_Muller(i*3+j,0.0,KbT/mass[i]);

  for (i=0;i<numatom;++i) 
    for (j=0;j<3;++j) 
      KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  return KE;
}

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

double MD_Propagetor_NH_MP1998_AAFF_Amber(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e,struct force *f) {
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

double MD_Propagetor_vV_NVE_AAFF_Amber(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f) {
  int i,j,k;
  double *frc,*frc_new,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frc_new=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_e[i*3+j]+(*f).f_LJ[i*3+j]+(*f).f_e_14[i*3+j]+(*f).f_LJ_14[i*3+j];

  ffL_calcffandforce(crd,numatom,e,f);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc_new[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_e[i*3+j]+(*f).f_LJ[i*3+j]+(*f).f_e_14[i*3+j]+(*f).f_LJ_14[i*3+j];

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=vel[i*3+j]+dt/mass[i]/2.0*(frc_new[i*3+j]+frc[i*3+j]);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]=crd[i*3+j]+dt*vel[i*3+j]+dt*dt/2.0/mass[i]*frc_new[i*3+j];

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  return KE;
}
