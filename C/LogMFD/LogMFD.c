
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <math.h>

#include "EF.h"
#include "FFL.h"
#include "PTL.h"

#include "TOPO.h"
#include "LA.h"

#include "MD_NHC_MP1996.h"
#include "MDrun.h"
#include "MD.h"

#include "Yacobian.h"

#include "LogMFD.h"

#define nys 3

double LogMFD_Propagetor_NH_Single_part_MP1996(double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int N,int nc,double wdt4[3],double wdt2[3]) {
  int i,j;

  double Scale;
  double AAAA;
  double G,KE;

  Scale=1.0;

  KE=0.0; for (i=0;i<N;++i) KE+=0.5*mass[i]*vel[i]*vel[i];
  G=(2.0*KE-NfKT)/Q;

  for (i=0;i<nc;++i) {
    for (j=0;j<nys;++j) {
      *V_zeta+=G*wdt4[j];
      AAAA=exp(-wdt2[j]*(*V_zeta));
      Scale=Scale*AAAA;
      G=(Scale*Scale*2.0*KE-NfKT)/Q;
      *zeta+=(*V_zeta)*wdt2[j];
      *V_zeta+=G*wdt4[j];
    }
  }

  for (i=0;i<N;++i) vel[i]=Scale*vel[i];

}

double LogMFD_Propagetor_NH_MP1998_1(double *crd, double *vel,double *mass,
				     double *zeta,double *V_zeta,double Q,
				     double NfKT,int N,double *KE,double *KEv,double *PEv,
				     double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
				     double *dF_dX, double F, double alpha, double gamma) {
  int i,j,k;
  double *frc;

  frc=(double *)gcemalloc(sizeof(double)*N);

  for (i=0;i<N;++i) {
    frc[i]=-(alpha*gamma/(alpha*F+1))*dF_dX[i];
  }

  LogMFD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,N,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<N;++i) vel[i]+=dt2*frc[i]/mass[i];
  for (i=0;i<N;++i) crd[i]+=dt*vel[i];

  *KE=0.0; for (i=0;i<N;++i) *KE+=0.5*mass[i]*vel[i]*vel[i];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return 0.0;
}

double LogMFD_Propagetor_NH_MP1998_2(double *crd, double *vel,double *mass,
				     double *zeta,double *V_zeta,double Q,
				     double NfKT,int N,double *KE,double *KEv,double *PEv,
				     double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
				     double *dF_dX, double F, double alpha, double gamma) {
  int i;
  double *frc;

  frc=(double *)gcemalloc(sizeof(double)*N);

  for (i=0;i<N;++i) {
    frc[i]=-(alpha*gamma/(alpha*F+1))*dF_dX[i];
  }

  for (i=0;i<N;++i) vel[i]+=dt2*frc[i]/mass[i];

  //////////////////////////////////////////////////////////////////////////////////////

  LogMFD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,N,nc,wdt4,wdt2);

  *KE=0.0; for (i=0;i<N;++i) *KE+=0.5*mass[i]*vel[i]*vel[i];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return 0.0;
}

double LogMFD_cF(double alpha, double gamma, double HMFD,
		 double NfKT, double KE,double KEv,double PEv) {
  double F,HMFD_KE_KEv_NfKT;

  HMFD_KE_KEv_NfKT=HMFD-KE-KEv-NfKT;

  F=(1.0/alpha)*(exp(HMFD_KE_KEv_NfKT/gamma)-1.0);

  return F;
}
