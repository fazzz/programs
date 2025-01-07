
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"

#include "FFL.h"

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

#include "TACCM_MDrun.h"
#include "TACCM_MD.h"
#include "MDrun.h"
#include "MD.h"
#include "MD_NHC_MP1996.h"

double runTACCM_MD_NHC_MP1998_Amber_AAFF(double *crd,double *vel, double *mass, int numatom,
					 double *zeta,double *V_zeta, double Q,
					 struct potential e, struct force f,
					 double T, double NfKT, int numstep, int interval,int *l,
					 double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
					 double *avePE, double *aveKE,double *aveT,
					 double *varPE, double *varKE,double *varT, double UNITT, double k_B,
					 struct my_netcdf_out_id_MCD nc_id_MCD,  FILE *outputfile,
					 //////////////// TACCM ///////////////////////
					 double *Z,double *velZ,double massZ,
					 double *zetaZ,double *V_zetaZ,
					 double TZ,double QZ,double NfKTZ,int numZ,
					 double Kapa,int **pairs,double pi,
					 double *avePEZ, double *aveKEZ,double *aveTZ,
					 double *varPEZ, double *varKEZ,double *varTZ, 
					 FILE *trjfileZ, FILE *trjfilTheta
					 //////////////// TACCM ///////////////////////
					 ) {
  int i,j,k;
  double PE=0.0,KE=0.0,Et,PEv,KEv;
  ///////////////// TACCM //////////////////////
  double *theta;
  double PEZ,KEZ,KEvZ,PEvZ,EtZ;
  double *frcZ;
  ///////////////// TACCM //////////////////////
  double summass,COM[3];
  double crd_nc[MAXATOM][3];

  *avePE=0.0; *aveKE=0.0; *aveT=0.0;
  *varPE=0.0; *varKE=0.0; *varT=0.0;
  *avePEZ=0.0; *aveKEZ=0.0; *aveTZ=0.0;
  *varPEZ=0.0; *varKEZ=0.0; *varTZ=0.0;

  ffL_calcffandforce(crd,numatom,&e,&f);

  frcZ=(double *)gcemalloc(sizeof(double)*numZ);
  theta=(double *)gcemalloc(sizeof(double)*numZ);
  TACCM_CTheta(crd,numatom,theta,numZ,pairs,pi);

  for (i=0;i<numstep;++i) {
    TACCM_CTheta(crd,numatom,theta,numZ,pairs,pi);

    //    TACCM_MD_Propagetor_NH_MP1998_AAFF_Amber(crd,vel,mass,&zeta,&V_zeta,Q,NfKT,numatom,&KE,&KEv,&PEv,dt,dt2,nc,wdt4,wdt2,&e,&f,Z,numZ,theta,Kapa,pairs,&PEZ,pi);
    TACCM_MD_Propagetor_NH_MP1998_AAFF_Amber(crd,vel,mass,zeta,V_zeta,Q,NfKT,numatom,&KE,&KEv,&PEv,dt,dt2,nc,wdt4,wdt2,&e,&f,Z,numZ,theta,Kapa,pairs,&PEZ,pi);

    //    TACCM_MD_Propagetor_NH_MP1998_Z(Z,velZ,massZ,theta,&zetaZ,&V_zetaZ,QZ,NfKTZ,numZ,&KEZ,&KEvZ,&PEvZ,dt,dt2,nc,wdt4,wdt2,Kapa,&PEZ,frcZ,pi);
    TACCM_MD_Propagetor_NH_MP1998_Z(Z,velZ,massZ,theta,zetaZ,V_zetaZ,QZ,NfKTZ,numZ,&KEZ,&KEvZ,&PEvZ,dt,dt2,nc,wdt4,wdt2,Kapa,&PEZ,frcZ,pi);
      
    if (i%interval==0) {
      KE=KE/UNITT;
      T=KE/((3*numatom)*k_B)*2.0;
      PEv=PEv/UNITT;
      KEv=KEv/UNITT;

      PE=0.5*e.p_e_t+0.5*e.p_LJ_t+0.5*e.p_e_14_t+0.5*e.p_LJ_14_t+e.p_d_t+e.p_a_t+e.p_b_t;
      Et=PE+KE+PEv+KEv;
      fprintf(outputfile,"%d %e %e %e %e %e %e %e ",i+1,PE,KE,KEv,PEv,Et,T);

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

      ///////////////// TACCM //////////////////////
      KEZ=KEZ/UNITT;
      TZ=KEZ/(numZ*k_B)*2.0;

      PEvZ=PEvZ/UNITT;
      KEvZ=KEvZ/UNITT;

      EtZ=PEZ+KEZ+PEvZ+KEvZ;
      fprintf(outputfile,"%d %e %e %e %e %e %e %e\n",i+1,PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);

      for (j=0;j<numZ;++j) fprintf(trjfileZ,"%e ",Z[j]);
      fprintf(trjfileZ,"\n");
      for (j=0;j<numZ;++j) fprintf(trjfilTheta,"%e ",theta[j]);
      fprintf(trjfilTheta,"\n");      
      ///////////////// TACCM //////////////////////
    }
  }

  return PEZ;
}

