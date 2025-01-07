
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"

#include "FFLc.h"

#include "MDrunb.h"
#include "MD.h"
#include "MDb_NHC_MP1996.h"

double runMDb_NHC_MP1998_Amber_AAFF(double *crd,double *vel, double *mass, int numatom,
				    double *zeta,double *V_zeta, double Q,
				    struct potential e, struct force f,
				    double T, double NfKT, int numstep, int interval,int *l,
				    double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
				    double *avePE, double *aveKE,double *aveT,
				    double *varPE, double *varKE,double *varT, double UNITT, double k_B,
				    struct my_netcdf_out_id_MCD nc_id_MCD,  FILE *outputfile,struct AmberParmL ap) {
  int i,j,k;
  double PE=0.0,KE=0.0,Et,PEv,KEv;
  double summass,COM[3];
  double crd_nc[MAXATOM][3];

  *avePE=0.0; *aveKE=0.0; *aveT=0.0;
  *varPE=0.0; *varKE=0.0; *varT=0.0;

  ffL_calcffandforce(crd,numatom,&e,&f);

  for (i=0;i<numstep;++i) {
    KE=MDb_Propagetor_NH_MP1998_AAFF_AMBER(crd,vel,mass,
					   zeta,V_zeta,Q,NfKT,numatom,&KEv,&PEv,dt,dt2,nc,wdt4,wdt2,&e,&f,ap);
      
    if (i%interval==0) {
      KE=KE/UNITT;
      T=KE/((3*numatom)*k_B)*2.0;
      PEv=PEv/UNITT;
      KEv=KEv/UNITT;

      PE=0.5*e.p_e_t+0.5*e.p_LJ_t+0.5*e.p_e_14_t+0.5*e.p_LJ_14_t+e.p_d_t+e.p_a_t+e.p_b_t;
      Et=PE+KE+PEv+KEv;
      fprintf(outputfile,"%d %e %e %e %e %e %e %e\n",i+1,PE,KE,KEv,PEv,Et,T);

      *avePE=(i*(*avePE)+PE)/(i+1);
      *varPE=(i*(*varPE)+PE*PE)/(i+1);
    
      *aveKE=(i*(*aveKE)+KE)/(i+1);
      *varKE=(i*(*varKE)+KE*KE)/(i+1);
    
      *aveT=(i*(*aveT)+T)/(i+1);
      *varT=(i*(*varT)+T*T)/(i+1);

      summass=0.0; for (j=0;j<numatom;++j) summass+=mass[j];
      for (j=0;j<3;++j) COM[j]=0.0;
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) COM[k]+=mass[j]*crd[j*3+k]/summass;
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) crd[j*3+k]-=COM[k];
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crd[j*3+k];
      myncL_put_crd_ene_MCD(nc_id_MCD,*l,crd_nc,e,0.0);
      ++(*l);
    }
  }

  return PE;
}
