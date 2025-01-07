
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MDb_NHC_MP1996.h"
#include "FFLc.h"
#include "PTLb.h"

#include "TOPO.h"
#include "LA.h"
#include "mymath.h"

#include "ABAb.h"
#include "TACCM.h"
#include "TACCM_MDb.h"
#include "FVDIHED.h"

#include "EF.h"

#include "RAND.h"
#include "BOXMULL.h"

#define UNITT 418.4070

#define ON  1
#define OFF 0

#define nys 3

double TACCMb_calc_eff_FF_MD(double *crd,int numatom, double *theta, double *Z,  int numZ,double Kapa, double **frcZ, int **pairs,double pi){
  int i,j,k;
  double PE=0.0;
  double *dvdpsi;
  double delta;

  int **pairs_temp;

  dvdpsi=(double *)gcemalloc(sizeof(double)*numZ);
  pairs_temp=(int **)gcemalloc(sizeof(int *)*numZ);
  for (i=0;i<numZ;++i) pairs_temp[i]=(int *)gcemalloc(sizeof(int)*4);

  for (i=0;i<numZ;++i) for (j=0;j<4;++j) pairs_temp[i][j]=pairs[i][j]-1;

  for (i=0;i<numZ;++i) {
    if ((delta=Z[i]-theta[i])>pi) delta-=2.0*pi;
    else if ((delta=Z[i]-theta[i])<-1.0*pi) delta+=2.0*pi;
    dvdpsi[i]=-Kapa*delta;
    PE+=0.5*Kapa*delta*delta;
  }

  FVDIHED_force_dihed(crd,numatom,frcZ,pairs_temp,dvdpsi,numZ);

  return PE;
}


double TACCMb_MD_Propagetor_NH_MP1998_AAFF_Amber(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KE,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f,double *Z,  int numZ,double *theta,double Kapa,int **pairs, double *PEZ,double pi, struct AmberParmL ap) {
  int i,j,k;
  double *frc;
  double **frcZ;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frcZ=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) frcZ[i]=(double *)gcemalloc(sizeof(double)*3);

  TACCM_CTheta(crd,numatom,theta,numZ,pairs,pi);
  *PEZ=TACCMb_calc_eff_FF_MD(crd,numatom,theta,Z,numZ,Kapa,frcZ,pairs,pi);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_e[i*3+j]+(*f).f_LJ[i*3+j]+(*f).f_e_14[i*3+j]+(*f).f_LJ_14[i*3+j]+frcZ[i][j];

  MDb_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];

  ffLc_calcffandforce(crd,numatom,e,f,ap);
  TACCM_CTheta(crd,numatom,theta,numZ,pairs,pi);
  *PEZ=TACCMb_calc_eff_FF_MD(crd,numatom,theta,Z,numZ,Kapa,frcZ,pairs,pi);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_e[i*3+j]+(*f).f_LJ[i*3+j]+(*f).f_e_14[i*3+j]+(*f).f_LJ_14[i*3+j]+frcZ[i][j];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MDb_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  *KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) *KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return 0.0;
}

double TACCMb_MD_Propagetor_NH_MP1998_Z(double *Z,double *velZ,double massZ,double *theta,
				       double *zeta,double *V_zeta,double Q,double NfKT,int numZ,
				       double *KE, double *KEv,double *PEv,
				       double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
				       double KZ,double *PEZ,double *f,double pi) {
  int i,j,k;

  TACCMb_MD_Propagetor_NH_Single_part_MP1996(velZ,massZ,zeta,V_zeta,Q,NfKT,numZ,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numZ;++i) velZ[i]+=dt2*f[i]/massZ;
  for (i=0;i<numZ;++i) {
    Z[i]+=dt*velZ[i];
    while (Z[i]>pi) Z[i]-=2.0*pi;
    while (Z[i]<=-pi) Z[i]+=2.0*pi;
  }
  *PEZ=TACCM_calc_eff_FF_Z(Z,numZ,theta,KZ,f,pi);
  for (i=0;i<numZ;++i) velZ[i]+=dt2*f[i]/massZ;
  //////////////////////////////////////////////////////////////////////////////////////

  TACCMb_MD_Propagetor_NH_Single_part_MP1996(velZ,massZ,zeta,V_zeta,Q,NfKT,numZ,nc,wdt4,wdt2);

  *KE=0.0; for (i=0;i<numZ;++i) *KE+=0.5*massZ*velZ[i]*velZ[i];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return 0.0;
}

double TACCMb_MD_Propagetor_vV_NVE_Z(double *Z, double *velZ,double massZ,double *theta,int numZ,
				    double dt, double KZ,double *PEZ,double *f,double pi) {
  int i,j,k;
  double *frc,*frc_new,KE;

  frc=(double *)gcemalloc(sizeof(double)*numZ);
  frc_new=(double *)gcemalloc(sizeof(double)*numZ);

  for (i=0;i<numZ;++i) frc[i]=f[i];

  *PEZ=TACCM_calc_eff_FF_Z(Z,numZ,theta,KZ,f,pi);
  for (i=0;i<numZ;++i) frc_new[i]=f[i];

  for (i=0;i<numZ;++i) velZ[i]=velZ[i]+dt/massZ/2.0*(frc_new[i]+frc[i]);
  for (i=0;i<numZ;++i) Z[i]=Z[i]+dt*velZ[i]+dt*dt/2.0/massZ*frc_new[i];

  KE=0.0; for (i=0;i<numZ;++i) KE+=0.5*massZ*velZ[i]*velZ[i];
  return KE;
}

double TACCMb_MD_Propagetor_NH_Single_part_MP1996(double *vel,double massZ,double *zeta,double *V_zeta,double Q,double NfKT,int numZ,int nc,double wdt4[3],double wdt2[3]) {
  int i,j;
  double Scale=1.0,AA,G,KE;

  KE=0.0; for (i=0;i<numZ;++i) KE+=0.5*massZ*vel[i]*vel[i];
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

  for (i=0;i<numZ;++i) vel[i]=Scale*vel[i];

}

double TACCMb_MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KE,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential_GOLMAA_PROTEINS2008 *e_GOLM,double *Z,  int numZ,double *theta,double Kapa,int **pairs, double *PEZ,double pi) {
  int i,j,k;
  double *frc;
  double **frcZ;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frcZ=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) frcZ[i]=(double *)gcemalloc(sizeof(double)*3);

  TACCM_CTheta(crd,numatom,theta,numZ,pairs,pi);
  *PEZ=TACCMb_calc_eff_FF_MD(crd,numatom,theta,Z,numZ,Kapa,frcZ,pairs,pi);
  GOLMAA_PROTEINS2008_ff_calcff_b(crd,numatom,e_GOLM);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_GOLM).f_b[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_d[i][j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j]+frcZ[i][j];

  MDb_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  GOLMAA_PROTEINS2008_ff_calcff_b(crd,numatom,e_GOLM);
  TACCM_CTheta(crd,numatom,theta,numZ,pairs,pi);
  *PEZ=TACCMb_calc_eff_FF_MD(crd,numatom,theta,Z,numZ,Kapa,frcZ,pairs,pi);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_GOLM).f_b[i][j]+(*e_GOLM).f_a[i][j]+(*e_GOLM).f_d[i][j]+(*e_GOLM).f_natatt[i][j]+(*e_GOLM).f_repul[i][j]+frcZ[i][j];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MDb_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  *KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) *KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return 0.0;
}
