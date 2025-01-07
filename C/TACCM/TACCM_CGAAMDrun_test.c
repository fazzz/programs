
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "FFL.h"

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

#include "TACCM_CGAAMDrun_test.h"
#include "TACCM_CGAAMDrun_test_CG.h"
#include "TACCM_MDrun.h"
#include "TACCM_MD.h"

#include "MD_NHC_MP1996.h"
#include "MDrun.h"
#include "MD.h"

double runTACCM_CGAA_MD_NHC_MP1998_test(// AA /////////////////////////////////////////////////////////
					double *crdAA,double *velAA, 
					double *zetaAA,double *V_zetaAA, double QAA,
					struct potential e, struct force f, double TAA, double NfKTAA,
					double *avePEAA, double *aveKEAA,double *aveTAA,
					double *varPEAA, double *varKEAA,double *varTAA,
					struct my_netcdf_out_id_MCD nc_id_MCDAA,  FILE *outputfileAA,
					// CG /////////////////////////////////////////////////////////
					double *crdCG,double *velCG, 
					double *zetaCG,double *V_zetaCG, double QCG,
					struct potential e_CG, struct force f_CG, double parameterCG,
					double TCG, double NfKTCG,
					double *avePECG, double *aveKECG,double *aveTCG,
					double *varPECG, double *varKECG,double *varTCG,
					struct my_netcdf_out_id_MCD nc_id_MCDCG,  FILE *outputfileCG,
					// Z  /////////////////////////////////////////////////////////
					double *Z,double *velZ,double massZ,
					double *zetaZ,double *V_zetaZ,
					double QZ,double NfKTZ,double TZ,
					int numZ,double KZAA, double KZCG,int **pairs,
					double *avePEZ, double *aveKEZ,double *aveTZ,
					double *varPEZ, double *varKEZ,double *varTZ, 
					FILE *trjfileZ, FILE *trjfilThetaAA, FILE *trjfilThetaCG,
					// CM  /////////////////////////////////////////////////////////
					double *mass, int numatom, int numstep, int interval,int *l,
					double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
					double UNITT, double k_B,double pi,
					double *PEZAA, double *PEZCG, double *PEZ) {
  int i,j,k;

  double PEAA=0.0,KEAA=0.0,EtAA,PEvAA,KEvAA;
  double PECG=0.0,KECG=0.0,EtCG,PEvCG,KEvCG;
  double KEZ,KEvZ,PEvZ,EtZ;

  double *thetaAA,*thetaCG,*frcZ,*fAA,*fCG;
  double summass,COM[3],crd_nc[MAXATOM][3];

  //  *avePEAA=0.0;
  //  *aveKEAA=0.0;
  //  *varPEAA=0.0;
  //  *varKEAA=0.0;
  //  *avePECG=0.0;
  //  *aveKECG=0.0;
  //  *varPECG=0.0;
  //  *varKECG=0.0;
  //  *avePEZ=0.0;
  //  *aveKEZ=0.0;
  //  *varPEZ=0.0;
  //  *varKEZ=0.0;

  //  *aveTAA=0.0;
  //  *varTAA=0.0;

  //  *aveTCG=0.0;
  //  *varTCG=0.0;

  //  *aveTZ=0.0;
  //  *varTZ=0.0;

  ffL_calcffandforce_AA(crdAA,numatom,&e,&f);
  ffL_calcffandforce_CG(crdCG,numatom,&e_CG,&f_CG,parameterCG);

  fAA=(double *)gcemalloc(sizeof(double)*numZ);
  fCG=(double *)gcemalloc(sizeof(double)*numZ);
  frcZ=(double *)gcemalloc(sizeof(double)*numZ);

  thetaAA=(double *)gcemalloc(sizeof(double)*numZ);
  thetaCG=(double *)gcemalloc(sizeof(double)*numZ);

  TACCM_CTheta(crdAA,numatom,thetaAA,numZ,pairs,pi);
  TACCM_CTheta(crdCG,numatom,thetaCG,numZ,pairs,pi);
  TACCM_calc_eff_FF_Z(Z,numZ,thetaAA,KZAA,fAA,pi);
  TACCM_calc_eff_FF_Z(Z,numZ,thetaCG,KZCG,fCG,pi);
  for (i=0;i<numZ;++i) frcZ[i]=fAA[i]+fCG[i];

  for (i=0;i<numstep;++i) {
    TACCM_CTheta(crdAA,numatom,thetaAA,numZ,pairs,pi);
    TACCM_CTheta(crdCG,numatom,thetaCG,numZ,pairs,pi);

    PEAA=TACCM_MD_Propagetor_NH_MP1998_AA_test(crdAA,velAA,mass,zetaAA,V_zetaAA,
					       QAA,NfKTAA,numatom,&KEAA,&KEvAA,&PEvAA,
					       dt,dt2,nc,wdt4,wdt2,
					       &e,&f,Z,numZ,thetaAA,KZAA,pairs,PEZAA,pi);

    PECG=TACCM_MD_Propagetor_NH_MP1998_CG_test(crdCG,velCG,mass,zetaCG,V_zetaCG,
					       QCG,NfKTCG,numatom,&KECG,&KEvCG,&PEvCG,
					       dt,dt2,nc,wdt4,wdt2,parameterCG,
					       &e_CG,&f_CG,Z,numZ,thetaCG,KZCG,pairs,PEZCG,pi);

    TACCM_CGAA_MD_Propagetor_NH_MP1998_Z_2(Z,velZ,massZ,thetaAA,thetaCG,zetaZ,V_zetaZ,
					   QZ,NfKTZ,numZ,&KEZ,&KEvZ,&PEvZ,
					   dt,dt2,nc,wdt4,wdt2,KZAA,KZCG,
					   PEZAA,PEZCG,PEZ,frcZ,pi);
      
    if (i%interval==0) {
      KEAA=KEAA/UNITT;     TAA=KEAA/((3*numatom)*k_B)*2.0;
      PEvAA=PEvAA/UNITT;   KEvAA=KEvAA/UNITT;

      KECG=KECG/UNITT;     TCG=KECG/((3*numatom)*k_B)*2.0;
      PEvCG=PEvCG/UNITT;  KEvCG=KEvCG/UNITT;

      PEAA=e.p_d_t+e.p_a_t+e.p_b_t;
      EtAA=PEAA+KEAA+PEvAA+KEvAA;
      fprintf(outputfileAA,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "
	                    ,i+1,PEAA,KEAA,KEvAA,PEvAA,EtAA,TAA);

      PECG=e_CG.p_d_t+e_CG.p_a_t+e_CG.p_b_t;
      EtCG=PECG+KECG+PEvCG+KEvCG;
      fprintf(outputfileCG,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "
	                  ,i+1,PECG,KECG,KEvCG,PEvCG,EtCG,TCG);

      //      *avePEAA=(i*(*avePEAA)+PEAA)/(i+1); *varPEAA=(i*(*varPEAA)+PEAA*PEAA)/(i+1);
      //      *aveKEAA=(i*(*aveKEAA)+KEAA)/(i+1); *varKEAA=(i*(*varKEAA)+KEAA*KEAA)/(i+1);
      //      *aveTAA=(i*(*aveTAA)+TAA)/(i+1);  *varTAA=(i*(*varTAA)+TAA*TAA)/(i+1);

      //      *avePECG=(i*(*avePECG)+PECG)/(i+1); *varPECG=(i*(*varPECG)+PECG*PECG)/(i+1);
      //      *aveKECG=(i*(*aveKECG)+KECG)/(i+1); *varKECG=(i*(*varKECG)+KECG*KECG)/(i+1);
      //      *aveTCG=(i*(*aveTCG)+TCG)/(i+1);  *varTCG=(i*(*varTCG)+TCG*TCG)/(i+1);

      summass=0.0; for (j=0;j<numatom;++j) summass+=mass[j];
      for (j=0;j<3;++j) COM[j]=0.0;
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) COM[k]+=mass[j]*crdAA[j*3+k]/summass;
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) crdAA[j*3+k]-=COM[k];
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crdAA[j*3+k];
      myncL_put_crd_ene_MCD(nc_id_MCDAA,*l,crd_nc,e,0.0);

      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crdCG[j*3+k];
      myncL_put_crd_ene_MCD(nc_id_MCDCG,*l,crd_nc,e,0.0);
      ++(*l);

      ///////////////// TACCM //////////////////////
      KEZ=KEZ/UNITT;      TZ=KEZ/(numZ*k_B)*2.0;

      PEvZ=PEvZ/UNITT;      KEvZ=KEvZ/UNITT;

      EtZ=*PEZ+KEZ+PEvZ+KEvZ;
      fprintf(outputfileAA,"%d %e %e %e %e %e %e %e\n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);

      for (j=0;j<numZ;++j) fprintf(trjfileZ,"%e ",Z[j]);
      fprintf(trjfileZ,"\n");
      for (j=0;j<numZ;++j) fprintf(trjfilThetaAA,"%e ",thetaAA[j]);
      fprintf(trjfilThetaAA,"\n");      
      for (j=0;j<numZ;++j) fprintf(trjfilThetaCG,"%e ",thetaCG[j]);
      fprintf(trjfilThetaCG,"\n");      
      ///////////////// TACCM //////////////////////
    }
  }

  return *PEZ;
}


/*********************************************************************************************************************/
/* double runTACCM_CGAA_MD_NHC_MP1998(// AA /////////////////////////////////////////////////////////		     */
/* 				   double *crdAA,double *velAA, 						     */
/* 				   double *zetaAA,double *V_zetaAA, double QAA,					     */
/* 				   struct potential e, struct force f, double TAA, double NfKTAA,		     */
/* 				   double *avePEAA, double *aveKEAA,double *aveTAA,				     */
/* 				   double *varPEAA, double *varKEAA,double *varTAA,				     */
/* 				   struct my_netcdf_out_id_MCD nc_id_MCDAA,  FILE *outputfileAA,		     */
/* 				   // CG /////////////////////////////////////////////////////////		     */
/* 				   double *crdCG,double *velCG, 						     */
/* 				   double *zetaCG,double *V_zetaCG, double QCG,					     */
/* 				   struct potential_GOLMAA_PROTEINS2008 e_GOLM, double TCG, double NfKTCG,	     */
/* 				   double *avePECG, double *aveKECG,double *aveTCG,				     */
/* 				   double *varPECG, double *varKECG,double *varTCG,				     */
/* 				   struct my_netcdf_out_id_MCD nc_id_MCDCG,  FILE *outputfileCG,		     */
/* 				   // Z  /////////////////////////////////////////////////////////		     */
/* 				   double *Z,double *velZ,double massZ,						     */
/* 				   double *zetaZ,double *V_zetaZ,						     */
/* 				   double QZ,double NfKTZ,double TZ,						     */
/* 				   int numZ,double KZAA, double KZCG,int **pairs,				     */
/* 				   double *avePEZ, double *aveKEZ,double *aveTZ,				     */
/* 				   double *varPEZ, double *varKEZ,double *varTZ, 				     */
/* 				   FILE *trjfileZ, FILE *trjfilThetaAA, FILE *trjfilThetaCG,			     */
/* 				   // CM  /////////////////////////////////////////////////////////		     */
/* 				   double *mass, int numatom, int numstep, int interval,int *l,			     */
/* 				   double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,			     */
/* 				   double UNITT, double k_B,double pi,						     */
/* 				   double *PEZAA, double *PEZCG, double *PEZ) {					     */
/*   int i,j,k;													     */
/* 														     */
/*   double PEAA=0.0,KEAA=0.0,EtAA,PEvAA,KEvAA;									     */
/*   double PECG=0.0,KECG=0.0,EtCG,PEvCG,KEvCG;									     */
/*   double KEZ,KEvZ,PEvZ,EtZ;											     */
/* 														     */
/*   double *thetaAA,*thetaCG,*frcZ,*fAA,*fCG;									     */
/*   double summass,COM[3],crd_nc[MAXATOM][3];									     */
/* 														     */
/*   *avePEAA=0.0; *aveKEAA=0.0; *aveTAA=0.0;									     */
/*   *varPEAA=0.0; *varKEAA=0.0; *varTAA=0.0;									     */
/*   *avePECG=0.0; *aveKECG=0.0; *aveTCG=0.0;									     */
/*   *varPECG=0.0; *varKECG=0.0; *varTCG=0.0;									     */
/*   *avePEZ=0.0; *aveKEZ=0.0; *aveTZ=0.0;									     */
/*   *varPEZ=0.0; *varKEZ=0.0; *varTZ=0.0;									     */
/* 														     */
/*   ffL_calcffandforce(crdAA,numatom,&e,&f);									     */
/* 														     */
/*   fAA=(double *)gcemalloc(sizeof(double)*numZ);								     */
/*   fCG=(double *)gcemalloc(sizeof(double)*numZ);								     */
/*   frcZ=(double *)gcemalloc(sizeof(double)*numZ);								     */
/* 														     */
/*   thetaAA=(double *)gcemalloc(sizeof(double)*numZ);								     */
/*   thetaCG=(double *)gcemalloc(sizeof(double)*numZ);								     */
/* 														     */
/*   TACCM_CTheta(crdAA,numatom,thetaAA,numZ,pairs,pi);								     */
/*   TACCM_CTheta(crdCG,numatom,thetaCG,numZ,pairs,pi);								     */
/*   TACCM_calc_eff_FF_Z(Z,numZ,thetaAA,KZAA,fAA,pi);								     */
/*   TACCM_calc_eff_FF_Z(Z,numZ,thetaCG,KZCG,fCG,pi);								     */
/*   for (i=0;i<numZ;++i) frcZ[i]=fAA[i]+fCG[i];								     */
/* 														     */
/*   for (i=0;i<numstep;++i) {											     */
/*     TACCM_CTheta(crdAA,numatom,thetaAA,numZ,pairs,pi);							     */
/*     TACCM_CTheta(crdCG,numatom,thetaCG,numZ,pairs,pi);							     */
/* 														     */
/*     PEAA=TACCM_MD_Propagetor_NH_MP1998_AAFF_Amber(crdAA,velAA,mass,zetaAA,V_zetaAA,				     */
/* 						  QAA,NfKTAA,numatom,&KEAA,&KEvAA,&PEvAA,			     */
/* 						  dt,dt2,nc,wdt4,wdt2,						     */
/* 						  &e,&f,Z,numZ,thetaAA,KZAA,pairs,PEZAA,pi);			     */
/* 														     */
/*     TACCM_MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008(crdCG,velCG,mass,zetaCG,V_zetaCG,			     */
/* 						      QCG,NfKTCG,numatom,&KECG,&KEvCG,&PEvCG,			     */
/* 						      dt,dt2,nc,wdt4,wdt2,					     */
/* 						      &e_GOLM,Z,numZ,thetaCG,KZCG,pairs,PEZCG,pi);		     */
/* 														     */
/*     TACCM_CGAA_MD_Propagetor_NH_MP1998_Z(Z,velZ,massZ,thetaAA,thetaCG,zetaZ,V_zetaZ,				     */
/* 					 QZ,NfKTZ,numZ,&KEZ,&KEvZ,&PEvZ,					     */
/* 					 dt,dt2,nc,wdt4,wdt2,KZAA,KZCG,PEZ,frcZ,pi);				     */
/*       													     */
/*     if (i%interval==0) {											     */
/*       KEAA=KEAA/UNITT;     TAA=KEAA/((3*numatom)*k_B)*2.0;							     */
/*       PEvAA=PEvAA/UNITT;   KEvAA=KEvAA/UNITT;								     */
/* 														     */
/*       KECG=KECG/UNITT;     TCG=KECG/((3*numatom)*k_B)*2.0;							     */
/*       PEvCG=PEvCG/UNITT;  KEvCG=KEvCG/UNITT;									     */
/* 														     */
/*       PEAA=0.5*e.p_e_t+0.5*e.p_LJ_t+0.5*e.p_e_14_t+0.5*e.p_LJ_14_t+e.p_d_t+e.p_a_t+e.p_b_t;			     */
/*       EtAA=PEAA+KEAA+PEvAA+KEvAA;										     */
/*       fprintf(outputfileAA,"%d %e %e %e %e %e %e %e ",i+1,PEAA,KEAA,KEvAA,PEvAA,EtAA,TAA);			     */
/* 														     */
/*       PECG=e_GOLM.p_t;											     */
/*       EtCG=PECG+KECG+PEvCG+KEvCG;										     */
/*       fprintf(outputfileCG,"%d %e %e %e %e %e %e %e ",i+1,PECG,KECG,KEvCG,PEvCG,EtCG,TCG);			     */
/* 														     */
/*       *avePEAA=(i*(*avePEAA)+PEAA)/(i+1); *varPEAA=(i*(*varPEAA)+PEAA*PEAA)/(i+1);				     */
/*       *aveKEAA=(i*(*aveKEAA)+KEAA)/(i+1); *varKEAA=(i*(*varKEAA)+KEAA*KEAA)/(i+1);				     */
/*       *aveTAA=(i*(*aveTAA)+TAA)/(i+1);  *varTAA=(i*(*varTAA)+TAA*TAA)/(i+1);					     */
/* 														     */
/*       *avePECG=(i*(*avePECG)+PECG)/(i+1); *varPECG=(i*(*varPECG)+PECG*PECG)/(i+1);				     */
/*       *aveKECG=(i*(*aveKECG)+KECG)/(i+1); *varKECG=(i*(*varKECG)+KECG*KECG)/(i+1);				     */
/*       *aveTCG=(i*(*aveTCG)+TCG)/(i+1);  *varTCG=(i*(*varTCG)+TCG*TCG)/(i+1);					     */
/* 														     */
/*       summass=0.0; for (j=0;j<numatom;++j) summass+=mass[j];							     */
/*       for (j=0;j<3;++j) COM[j]=0.0;										     */
/*       for (j=0;j<numatom;++j)  for (k=0;k<3;++k) COM[k]+=mass[j]*crdAA[j*3+k]/summass;			     */
/*       for (j=0;j<numatom;++j)  for (k=0;k<3;++k) crdAA[j*3+k]-=COM[k];					     */
/*       for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crdAA[j*3+k];					     */
/*       myncL_put_crd_ene_MCD(nc_id_MCDAA,*l,crd_nc,e,0.0);							     */
/* 														     */
/*       for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crdCG[j*3+k];					     */
/*       myncL_put_crd_ene_MCD(nc_id_MCDCG,*l,crd_nc,e,0.0);							     */
/*       ++(*l);												     */
/* 														     */
/*       ///////////////// TACCM //////////////////////								     */
/*       KEZ=KEZ/UNITT;      TZ=KEZ/(numZ*k_B)*2.0;								     */
/* 														     */
/*       PEvZ=PEvZ/UNITT;      KEvZ=KEvZ/UNITT;									     */
/* 														     */
/*       EtZ=*PEZ+KEZ+PEvZ+KEvZ;										     */
/*       fprintf(outputfileAA,"%d %e %e %e %e %e %e %e\n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);			     */
/* 														     */
/*       for (j=0;j<numZ;++j) fprintf(trjfileZ,"%e ",Z[j]);							     */
/*       fprintf(trjfileZ,"\n");										     */
/*       for (j=0;j<numZ;++j) fprintf(trjfilThetaAA,"%e ",thetaAA[j]);						     */
/*       fprintf(trjfilThetaAA,"\n");      									     */
/*       for (j=0;j<numZ;++j) fprintf(trjfilThetaCG,"%e ",thetaCG[j]);						     */
/*       fprintf(trjfilThetaCG,"\n");      									     */
/*       ///////////////// TACCM //////////////////////								     */
/*     }													     */
/*   }														     */
/* 														     */
/*   return *PEZ;												     */
/* }														     */
/*********************************************************************************************************************/


double runMD_NHC_MP1998_Amber_AAFF_test(double *crd,double *vel, double *mass, int numatom,
					double *zeta,double *V_zeta, double Q,
					struct potential e, struct force f,
					double T, double NfKT, int numstep, int interval,int *l,
					double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
					double *avePE, double *aveKE,double *aveT,
					double *varPE, double *varKE,double *varT, double UNITT, double k_B,
					struct my_netcdf_out_id_MCD nc_id_MCD,  FILE *outputfile) {
  int i,j,k;
  double PE=0.0,KE=0.0,Et,PEv,KEv;
  double summass,COM[3];
  double crd_nc[MAXATOM][3];

  *avePE=0.0; *aveKE=0.0; *aveT=0.0;
  *varPE=0.0; *varKE=0.0; *varT=0.0;

  ffL_calcffandforce(crd,numatom,&e,&f);

  for (i=0;i<numstep;++i) {
    KE=MD_Propagetor_NH_MP1998_AAFF_Amber_test(crd,vel,mass,
					       zeta,V_zeta,Q,NfKT,numatom,&KEv,&PEv,dt,dt2,nc,wdt4,wdt2,&e,&f);
      
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

double runTACCM_MD_NHC_MP1998_Amber_AAFF_test(double *crd,double *vel, double *mass, int numatom,
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
    TACCM_MD_Propagetor_NH_MP1998_AA_test(crd,vel,mass,
					  zeta,V_zeta,Q,NfKT,numatom,
					  &KE,&KEv,&PEv,
					  dt,dt2,nc,wdt4,wdt2,&e,&f,Z,numZ,theta,Kapa,pairs,&PEZ,pi);

    //    TACCM_MD_Propagetor_NH_MP1998_Z(Z,velZ,massZ,theta,&zetaZ,&V_zetaZ,QZ,NfKTZ,numZ,&KEZ,&KEvZ,&PEvZ,dt,dt2,nc,wdt4,wdt2,Kapa,&PEZ,frcZ,pi);
    TACCM_MD_Propagetor_NH_MP1998_Z(Z,velZ,massZ,theta,
				    zetaZ,V_zetaZ,QZ,NfKTZ,numZ,
				    &KEZ,&KEvZ,&PEvZ,
				    dt,dt2,nc,wdt4,wdt2,Kapa,&PEZ,frcZ,pi);
      
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

double CE_TACCM_CGAA_testb(double *crdAA,double *crdCG,double *Z, int numatom,int numZ,
			   double KZAA,double KZCG,int **pairs,double pi,
			   double *EAA,double *ECG,double *EZ){
  int i,j;
  double *thetaAA,*thetaCG;
  double delta;

  thetaAA=(double *)gcemalloc(sizeof(double)*numZ);
  thetaCG=(double *)gcemalloc(sizeof(double)*numZ);

  TACCM_CTheta(crdAA,numatom,thetaAA,numZ,pairs,pi);
  TACCM_CTheta(crdCG,numatom,thetaCG,numZ,pairs,pi);

  *EAA=0.0;
  for (i=0;i<numZ;++i) {
    if ((delta=Z[i]-thetaAA[i])>pi) delta-=2.0*pi;
    else if ((delta=Z[i]-thetaAA[i])<-1.0*pi) delta+=2.0*pi;
    *EAA+=0.5*KZAA*delta*delta;
  }

  *ECG=0.0;
  for (i=0;i<numZ;++i) {
    if ((delta=Z[i]-thetaCG[i])>pi) delta-=2.0*pi;
    else if ((delta=Z[i]-thetaCG[i])<-1.0*pi) delta+=2.0*pi;
    *ECG+=0.5*KZCG*delta*delta;
  }

  *EZ=(*EAA)+(*ECG);
}

double MD_Propagetor_NH_MP1998_AAFF_Amber_test(double *crd,double *vel,double *mass,double *zeta,double *V_zeta,double Q,double NfKT,int numatom,double *KEv,double *PEv,double dt,double dt2,int nc,double wdt4[3],double wdt2[3], struct potential *e, struct force *f) {
  int i,j,k;
  double *frc,KE;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j];
  
  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];
  ffL_calcffandforce_AA(crd,numatom,e,f);
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return KE;
}

double TACCM_MD_Propagetor_NH_MP1998_AA_test(double *crd,double *vel,double *mass,
					     double *zeta,double *V_zeta,double Q,
					     double NfKT,int numatom,double *KE,double *KEv,double *PEv,
					     double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
					     struct potential *e, struct force *f,
					     double *Z,  int numZ,double *theta,double Kapa,
					     int **pairs, double *PEZ,double pi) {
  int i,j,k;
  double *frc;
  double **frcZ;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frcZ=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) frcZ[i]=(double *)gcemalloc(sizeof(double)*3);

  TACCM_CTheta(crd,numatom,theta,numZ,pairs,pi);
  *PEZ=TACCM_calc_eff_FF_MD(crd,numatom,theta,Z,numZ,Kapa,frcZ,pairs,pi);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+frcZ[i][j];

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];

  ffL_calcffandforce_AA(crd,numatom,e,f);
  TACCM_CTheta(crd,numatom,theta,numZ,pairs,pi);
  *PEZ=TACCM_calc_eff_FF_MD(crd,numatom,theta,Z,numZ,Kapa,frcZ,pairs,pi);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+frcZ[i][j];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  *KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) *KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return 0.0;
}
