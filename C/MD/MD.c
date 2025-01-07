
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "RAND.h"
#include "BOXMULL.h"
#include "MD.h"

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

double MD_Propagetor_Iso(double *crd,double *vel,double *mass,int numatom,double IsoCoff,double dt,double *KE,double *PE) {
  int i,j,k;
  double *velp,*crd_h,*frc;
  double sum;
  struct potential e;
  struct force f;

  velp=(double *)gcemalloc(sizeof(double)*numatom*3);
  crd_h=(double *)gcemalloc(sizeof(double)*numatom*3);
  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j)
      crd_h[i*3+j]=crd[i*3+j]+dt*vel[i*3+j];

  ffL_set_calcffandforce(&e,&f);
  for (i=0;i<numatom*3;++i) {
    f.f_b[i]=0.0;    f.f_a[i]=0.0;    f.f_d[i]=0.0;
    f.f_e[i]=0.0;    f.f_LJ[i]=0.0;   f.f_e_14[i]=0.0;   f.f_LJ_14[i]=0.0;
  }
  ffL_calcffandforce(crd_h,numatom,&e,&f);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-f.f_b[i*3+j]+f.f_a[i*3+j]+f.f_d[i*3+j]+f.f_e[i*3+j]+f.f_LJ[i*3+j]+f.f_e_14[i*3+j]+f.f_LJ_14[i*3+j];
  *(PE)=e.p_b_t+e.p_a_t+e.p_d_t+0.5*e.p_e_t+0.5*e.p_LJ_t+0.5*e.p_e_14_t+0.5*e.p_LJ_14_t;

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j)
      velp[i*3+j]=vel[i*3+j]+0.5*dt/mass[i]*frc[i*3+j];

  sum=0.0;
  for (i=0;i<numatom;++i)for (j=0;j<3;++j) sum+=velp[i*3+j]*velp[i*3+j];
  sum=sqrt(sum);

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      vel[i*3+j]=IsoCoff*velp[i*3+j]/sum;
      crd[i*3+j]+=dt*vel[i*3+j];
    }
  }

  *(KE)=0.0;
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) *(KE)+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
}

double MD_Propagetor_Iso_JCP2003(double *crd,double *vel,double *mass,int numatom,double K,double dt,struct potential *e, struct force *f) {
  int i,j,k;
  double *frc,KE;
  double a,b,s,ds_dt;
  double UNITT=418.4070;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_e[i*3+j]+(*f).f_LJ[i*3+j]+(*f).f_e_14[i*3+j]+(*f).f_LJ_14[i*3+j];

  a=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) a+=(frc[i*3+j]*vel/*0*/[i*3+j])/2.0/K;
  b=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) b+=(frc[i*3+j]*frc[i*3+j]/mass[i])/2.0/K;
  s=a/b*(cosh(0.5*dt*sqrt(b))-1.0)+1.0/sqrt(b)*(sinh(0.5*dt*sqrt(b)));
  ds_dt=a/sqrt(b)*sinh(0.5*dt*sqrt(b))+cosh(0.5*dt*sqrt(b));

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=(mass[i]*vel[i*3+j]+frc[i*3+j]*s)/ds_dt/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j]/*/mass[i]*/;

  ffL_calcffandforce(crd,numatom,e,f);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_e[i*3+j]+(*f).f_LJ[i*3+j]+(*f).f_e_14[i*3+j]+(*f).f_LJ_14[i*3+j];

  a=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) a+=(frc[i*3+j]*vel/*0*/[i*3+j])/2.0/K;
  b=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) b+=(frc[i*3+j]*frc[i*3+j]/mass[i])/2.0/K;
  s=a/b*(cosh(0.5*dt*sqrt(b))-1.0)+1.0/sqrt(b)*(sinh(0.5*dt*sqrt(b)));
  ds_dt=a/sqrt(b)*sinh(0.5*dt*sqrt(b))+cosh(0.5*dt*sqrt(b));

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=(vel[i*3+j]*mass[i]+frc[i*3+j]*s)/ds_dt/mass[i];

  KE=0.0;
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  return KE;
}

double MD_Propagetor_Iso_JCP2003_GOLMAA_JCTC2011(double *crd,double *vel,double *mass,int numatom,double K,double dt,struct potential *e, struct force *f, struct potential_GOLMAA_hybrid *e_GOLM) {
  int i,j,k;
  double *frc,KE;
  double a,b,s,ds_dt;
  double UNITT=418.4070;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];

  a=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) a+=(frc[i*3+j]*vel/*0*/[i*3+j])/2.0/K;
  b=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) b+=(frc[i*3+j]*frc[i*3+j]/mass[i])/2.0/K;
  s=a/b*(cosh(0.5*dt*sqrt(b))-1.0)+1.0/sqrt(b)*(sinh(0.5*dt*sqrt(b)));
  ds_dt=a/sqrt(b)*sinh(0.5*dt*sqrt(b))+cosh(0.5*dt*sqrt(b));

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=(mass[i]*vel[i*3+j]+frc[i*3+j]*s)/ds_dt/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j]/*/mass[i]*/;

  ffL_calcffandforce_14DAB_woH(crd,numatom,e,f);
  //  ffL_calcffandforce_14DAB_woH(crd,numatom,e,f);
  GOLMAA_hyb_ff_calcff(crd,numatom,e_GOLM);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];

  a=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) a+=(frc[i*3+j]*vel/*0*/[i*3+j])/2.0/K;
  b=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) b+=(frc[i*3+j]*frc[i*3+j]/mass[i])/2.0/K;
  s=a/b*(cosh(0.5*dt*sqrt(b))-1.0)+1.0/sqrt(b)*(sinh(0.5*dt*sqrt(b)));
  ds_dt=a/sqrt(b)*sinh(0.5*dt*sqrt(b))+cosh(0.5*dt*sqrt(b));

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=(vel[i*3+j]*mass[i]+frc[i*3+j]*s)/ds_dt/mass[i];

  KE=0.0;
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  return KE;
}

double MD_Propagetor_vV_NVE_GOLMAA_JCTC2011(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f, struct potential_GOLMAA_hybrid *e_GOLM) {
  int i,j,k;
  double *frc,*frc_new,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frc_new=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];

  ffL_calcffandforce_14vdWDAB_woH(crd,numatom,e,f);
  GOLMAA_hyb_ff_calcff(crd,numatom,e_GOLM);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc_new[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=vel[i*3+j]+dt/mass[i]/2.0*(frc_new[i*3+j]+frc[i*3+j]);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]=crd[i*3+j]+dt*vel[i*3+j]+dt*dt/2.0/mass[i]*frc_new[i*3+j];

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  return KE;
}

double MD_Propagetor_vV_NVE_14LJdab(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f) {
  int i,j,k;
  double *frc,*frc_new,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frc_new=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j];

  ffL_calcffandforce_14vdWDAB_woH(crd,numatom,e,f);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc_new[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j];

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=vel[i*3+j]+dt/mass[i]/2.0*(frc_new[i*3+j]+frc[i*3+j]);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]=crd[i*3+j]+dt*vel[i*3+j]+dt*dt/2.0/mass[i]*frc_new[i*3+j];

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  return KE;
}

double MD_Propagetor_vV_NVE_GOLMAA_PROTEINS2008(double *crd,double *vel,double *mass,int numatom,double dt,struct potential_GOLMAA_PROTEINS2008 *e_GOLM) {
  int i,j,k;
  double *frc,*frc_new,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frc_new=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_GOLM).f_b[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_d[i][j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];

  GOLMAA_PROTEINS2008_ff_calcff(crd,numatom,e_GOLM);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc_new[i*3+j]=(*e_GOLM).f_b[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_d[i][j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=vel[i*3+j]+dt/mass[i]/2.0*(frc_new[i*3+j]+frc[i*3+j]);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]=crd[i*3+j]+dt*vel[i*3+j]+dt*dt/2.0/mass[i]*frc_new[i*3+j];

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  return KE;
}

double MD_Propagetor_vV_NVE_GOLMAA_PROTEINS2008_b(double *crd,double *vel,double *mass,int numatom,double dt,struct potential_GOLMAA_PROTEINS2008 *e_GOLM) {
  int i,j,k;
  double *frc,*frc_new,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frc_new=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_GOLM).f_b[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_d[i][j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];

  GOLMAA_PROTEINS2008_ff_calcff_b(crd,numatom,e_GOLM);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc_new[i*3+j]=(*e_GOLM).f_b[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_d[i][j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j];

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=vel[i*3+j]+dt/mass[i]/2.0*(frc_new[i*3+j]+frc[i*3+j]);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]=crd[i*3+j]+dt*vel[i*3+j]+dt*dt/2.0/mass[i]*frc_new[i*3+j];

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  return KE;
}

double MD_Propagetor_vV_NVE_GOLMAA_MB_PROTEINS2008(double *crd,double *vel,double *mass,int numatom,double dt,double de, double d2,struct potential_GOLMAA_MB_PROTEINS2008 *e_GOLM) {
  int i,j,k;
  double *frc,*frc_new,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frc_new=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) frc[i*3+j]=(*e_GOLM).f_MB[i][j];

  GOLMAA_MB_PROTEINS2008_ff_calcff(crd,numatom,de,d2,e_GOLM);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc_new[i*3+j]=(*e_GOLM).f_MB[i][j];

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=vel[i*3+j]+dt/mass[i]/2.0*(frc_new[i*3+j]+frc[i*3+j]);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]=crd[i*3+j]+dt*vel[i*3+j]+dt*dt/2.0/mass[i]*frc_new[i*3+j];

  KE=0.0; for (i=0;i<numatom;++i) 
	    for (j=0;j<3;++j) 
	      KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  return KE;
}

double MD_Propagetor_vV_NVE_GOLMAA_JCTC2011_wonat(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f, struct potential_GOLMAA_hybrid *e_GOLM) {
  int i,j,k;
  double *frc,*frc_new,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frc_new=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j]+(*e_GOLM).f_repul[i][j];

  ffL_calcffandforce_14vdWDAB_woH(crd,numatom,e,f);
  GOLMAA_hyb_ff_calcff(crd,numatom,e_GOLM);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc_new[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j]+(*e_GOLM).f_repul[i][j];

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=vel[i*3+j]+dt/mass[i]/2.0*(frc_new[i*3+j]+frc[i*3+j]);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]=crd[i*3+j]+dt*vel[i*3+j]+dt*dt/2.0/mass[i]*frc_new[i*3+j];

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  return KE;
}

double MD_Propagetor_vV_NVE_14LJdba(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f) {
  int i,j,k;
  double *frc,*frc_new,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frc_new=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j];

  ffL_calcffandforce_14vdWDAB_woH(crd,numatom,e,f);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc_new[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j];

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=vel[i*3+j]+dt/mass[i]/2.0*(frc_new[i*3+j]+frc[i*3+j]);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]=crd[i*3+j]+dt*vel[i*3+j]+dt*dt/2.0/mass[i]*frc_new[i*3+j];

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  return KE;
}

double MD_Propagetor_Iso_JCP2003_GOLM_Clementi(double *crd,double *vel,double *mass,int numatom,double K,double dt, struct potential_GOLM_Clementi *e_GOLM) {
  int i,j,k;
  double *frc,KE;
  double a,b,s,ds_dt;
  double UNITT=418.4070;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_GOLM).f_t[i][j];

  a=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) a+=(frc[i*3+j]*vel/*0*/[i*3+j])/2.0/K;
  b=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) b+=(frc[i*3+j]*frc[i*3+j]/mass[i])/2.0/K;
  s=a/b*(cosh(0.5*dt*sqrt(b))-1.0)+1.0/sqrt(b)*(sinh(0.5*dt*sqrt(b)));
  ds_dt=a/sqrt(b)*sinh(0.5*dt*sqrt(b))+cosh(0.5*dt*sqrt(b));

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=(mass[i]*vel[i*3+j]+frc[i*3+j]*s)/ds_dt/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j]/*/mass[i]*/;

  GOLM_Clementi_ff_calcff(crd,numatom,e_GOLM);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_GOLM).f_t[i][j];

  a=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) a+=(frc[i*3+j]*vel/*0*/[i*3+j])/2.0/K;
  b=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) b+=(frc[i*3+j]*frc[i*3+j]/mass[i])/2.0/K;
  s=a/b*(cosh(0.5*dt*sqrt(b))-1.0)+1.0/sqrt(b)*(sinh(0.5*dt*sqrt(b)));
  ds_dt=a/sqrt(b)*sinh(0.5*dt*sqrt(b))+cosh(0.5*dt*sqrt(b));

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=(vel[i*3+j]*mass[i]+frc[i*3+j]*s)/ds_dt/mass[i];

  KE=0.0;
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  return KE;
}

double MD_Propagetor_vV_NVE_GOLM_Clementi(double *crd,double *vel,double *mass,int numatom,double dt,struct potential_GOLM_Clementi *e_GOLM) {
  int i,j,k;
  double *frc,*frc_new,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frc_new=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_GOLM).f_t[i][j];;

  GOLM_Clementi_ff_calcff(crd,numatom,e_GOLM);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc_new[i*3+j]=(*e_GOLM).f_t[i][j];;

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=vel[i*3+j]+dt/mass[i]/2.0*(frc_new[i*3+j]+frc[i*3+j]);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]=crd[i*3+j]+dt*vel[i*3+j]+dt*dt/2.0/mass[i]*frc_new[i*3+j];

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  return KE;
}

double MD_Propagetor_vV_NVE(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f) {
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

double MD_Propagetor_vV_NVE_wc(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f,int bflag,int aflag,int dflag,int eflag, int LJflag,int e14flag,int LJ14flag ) {
  int i,j,k;
  double *frc,*frc_new,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frc_new=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      frc[i*3+j]=0.0;
      if(bflag==ON)   frc[i*3+j]-=(*f).f_b[i*3+j];
      if(aflag==ON)  frc[i*3+j]+=(*f).f_a[i*3+j];
      if(dflag==ON)  frc[i*3+j]+=(*f).f_d[i*3+j];
      if (e14flag == ON)  frc[i*3+j]+=(*f).f_e_14[i*3+j];
      if (LJ14flag== ON)  frc[i*3+j]+=(*f).f_LJ_14[i*3+j];
      if (eflag == ON)  frc[i*3+j]+=(*f).f_e[i*3+j];
      if (LJflag== ON)  frc[i*3+j]+=(*f).f_LJ[i*3+j];
    }
  }

  ffL_calcffandforce(crd,numatom,e,f);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      frc_new[i*3+j]=0.0;
      if(bflag==ON)  frc_new[i*3+j]-=(*f).f_b[i*3+j];
      if(aflag==ON)  frc_new[i*3+j]+=(*f).f_a[i*3+j];
      if(dflag==ON)  frc_new[i*3+j]+=(*f).f_d[i*3+j];
      if (e14flag == ON)  frc_new[i*3+j]+=(*f).f_e_14[i*3+j];
      if (LJ14flag== ON)  frc_new[i*3+j]+=(*f).f_LJ_14[i*3+j];
      if (eflag == ON)  frc_new[i*3+j]+=(*f).f_e[i*3+j];
      if (LJflag== ON)  frc_new[i*3+j]+=(*f).f_LJ[i*3+j];
    }
  }

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=vel[i*3+j]+dt/mass[i]/2.0*(frc_new[i*3+j]+frc[i*3+j]);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]=crd[i*3+j]+dt*vel[i*3+j]+dt*dt/2.0/mass[i]*frc_new[i*3+j];

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  return KE;
}

double MD_Propagetor_vV_NVE_FASYS(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f) {
  int i,j,k;
  double *frc,*frc_new,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frc_new=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j];

  ffL_calcffandforce(crd,numatom,e,f);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc_new[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j];

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=vel[i*3+j]+dt/mass[i]/2.0*(frc_new[i*3+j]-frc[i*3+j]);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]=crd[i*3+j]+dt*vel[i*3+j]+dt*dt/2.0/mass[i]*frc_new[i*3+j];

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  return KE;
}

double MD_Propagetor_vV_NVE_wflag(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f,int bodflag,int angflag,int dihflag,int LJ14flag,int es14flag,int LJflag,int esflag) {
  int i,j,k;
  double *frc,*frc_new,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frc_new=(double *)gcemalloc(sizeof(double)*numatom*3);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      frc[i*3+j]=0.0;
      frc_new[i*3+j]=0.0;
    }
  }

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

  ffL_calcffandforce(crd,numatom,e,f);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      if (esflag==ON)   frc_new[i*3+j]+=(*f).f_e[i*3+j];
      if (LJflag==ON)   frc_new[i*3+j]+=(*f).f_LJ[i*3+j];
      if (es14flag==ON) frc_new[i*3+j]+=(*f).f_e_14[i*3+j];
      if (LJ14flag==ON) frc_new[i*3+j]+=(*f).f_LJ_14[i*3+j];
      if (dihflag==ON)  frc_new[i*3+j]+=(*f).f_d[i*3+j];
      if (angflag==ON)  frc_new[i*3+j]+=(*f).f_a[i*3+j];
      if (bodflag==ON)  frc_new[i*3+j]-=(*f).f_b[i*3+j];
    }
  }

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=vel[i*3+j]+dt/mass[i]/2.0*(frc_new[i*3+j]+frc[i*3+j]);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]=crd[i*3+j]+dt*vel[i*3+j]+dt*dt/2.0/mass[i]*frc_new[i*3+j];

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  return KE;
}

double MD_Propagetor_vV_NVE_woH_wflag(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f,int bodflag,int angflag,int dihflag,int LJ14flag,int es14flag,int LJflag,int esflag) {
  int i,j,k;
  double *frc,*frc_new,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frc_new=(double *)gcemalloc(sizeof(double)*numatom*3);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      frc[i*3+j]=0.0;
      frc_new[i*3+j]=0.0;
    }
  }

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

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=vel[i*3+j]+dt/mass[i]/2.0*(frc_new[i*3+j]+frc[i*3+j]);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]=crd[i*3+j]+dt*vel[i*3+j]+dt*dt/2.0/mass[i]*frc_new[i*3+j];

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  return KE;
}

double MD_Propagetor_vV_NVE_AAFF_Amber(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f) {
  int i,j,k;
  double *frc,*frc_new,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frc_new=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_e_14[i*3+j]+(*f).f_LJ_14[i*3+j]+(*f).f_e[i*3+j]+(*f).f_LJ[i*3+j];

  ffL_calcffandforce(crd,numatom,e,f);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc_new[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_e_14[i*3+j]+(*f).f_LJ_14[i*3+j]+(*f).f_e[i*3+j]+(*f).f_LJ[i*3+j];

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=vel[i*3+j]+dt/mass[i]/2.0*(frc_new[i*3+j]+frc[i*3+j]);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]=crd[i*3+j]+dt*vel[i*3+j]+dt*dt/2.0/mass[i]*frc_new[i*3+j];

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  return KE;
}

double MD_Propagetor_vV_NVE_UMB(double *crd,double *vel,double *mass,int numatom,double dt,struct potential *e, struct force *f,int *pairp,int nump,double *fcp,double *dih_equ,double *pUMB,double **fUMB) {
  int i,j,k;
  double *frc,*frc_new,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frc_new=(double *)gcemalloc(sizeof(double)*numatom*3);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      frc[i*3+j]=0.0;
      frc_new[i*3+j]=0.0;
    }
  }

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j]+(*f).f_e_14[i*3+j]+(*f).f_LJ[i*3+j]+(*f).f_e[i*3+j]+fUMB[i][j];
    }
  }

  ffL_calcffandforce(crd,numatom,e,f);
  UMB_calc_dihetype_ff(crd,numatom,pairp,nump,fcp,dih_equ,pUMB,fUMB);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      frc_new[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_LJ_14[i*3+j]+(*f).f_e_14[i*3+j]+(*f).f_LJ[i*3+j]+(*f).f_e[i*3+j]+fUMB[i][j];
    }
  }

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]=vel[i*3+j]+dt/mass[i]/2.0*(frc_new[i*3+j]+frc[i*3+j]);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]=crd[i*3+j]+dt*vel[i*3+j]+dt*dt/2.0/mass[i]*frc_new[i*3+j];

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  return KE;
}
