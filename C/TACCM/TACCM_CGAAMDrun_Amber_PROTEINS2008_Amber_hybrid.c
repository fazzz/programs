
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "FFL.h"

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

#include "TACCM_CGAAMDrun_Amber_PROTEINS2008_Amber_hybrid.h"
//#include "TACCM_CGAAMDrun_test_CG.h"
#include "TACCM_MDrun.h"
#include "TACCM_MD.h"

#include "MD_NHC_MP1996.h"
#include "MDrun.h"
#include "MD.h"

double runTACCM_CGAA_MD_NHC_MP1998_Amber_PROTEINS2008_Amber_hybrid(// AA ////////////////////////////////////////////
								   double *crdAA,double *velAA, 
								   double *zetaAA,double *V_zetaAA, double QAA,
								   struct potential e, struct force f,
								   double TAA, double NfKTAA,
								   double *avePEAA, double *aveKEAA,double *aveTAA,
								   double *varPEAA, double *varKEAA,double *varTAA,
								   struct my_netcdf_out_id_MCD nc_id_MCDAA,
								   FILE *outputfileAA,
								   // CG ////////////////////////////////////////////
								   double *crdCG,double *velCG, 
								   double *zetaCG,double *V_zetaCG, double QCG,
								   struct potential_GOLMAA_PROTEINS2008 e_CG, 
								   double TCG, double NfKTCG,
								   double *avePECG, double *aveKECG,double *aveTCG,
								   double *varPECG, double *varKECG,double *varTCG,
								   struct my_netcdf_out_id_MCD nc_id_MCDCG,
								   FILE *outputfileCG,
								   // Z  ////////////////////////////////////////////
								   double *Z,double *velZ,double massZ,
								   double *zetaZ,double *V_zetaZ,
								   double QZ,double NfKTZ,double TZ,
								   int numZ,double KZAA, double KZCG,int **pairs,
								   double *avePEZ, double *aveKEZ,double *aveTZ,
								   double *varPEZ, double *varKEZ,double *varTZ, 
								   FILE *trjfileZ, FILE *trjfilThetaAA, 
								   FILE *trjfilThetaCG,
								   // CM  //////////////////////////////////////////
								   double *mass, int numatom, int numheavyatom,
								   int numstep, int interval,int *l,
								   double dt,double dt2,
								   double wdt2[3] ,double wdt4[3] ,int nc,
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

  ffL_calcffandforce(crdAA,numatom,&e,&f);
  //  ffL_calcffandforce_CG(crdCG,numatom,&e_CG,&f_CG,parameterCG);
  GOLMAA_PROTEINS2008_Amber_hybrid_ff_calcff_b(crdCG,numatom,&e_CG);

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

    PEAA=TACCM_MD_Propagetor_NH_MP1998_AAFF_Amber(crdAA,velAA,mass,zetaAA,V_zetaAA,
						  QAA,NfKTAA,numatom,&KEAA,&KEvAA,&PEvAA,
						  dt,dt2,nc,wdt4,wdt2,
						  &e,&f,Z,numZ,thetaAA,KZAA,pairs,PEZAA,pi);

    PECG=TACCM_MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008_Amber_hybrid(crdCG,velCG,mass,zetaCG,V_zetaCG,
									QCG,NfKTCG,numatom,&KECG,&KEvCG,&PEvCG,
									dt,dt2,nc,wdt4,wdt2,
									&e_CG,Z,numZ,thetaCG,KZCG,pairs,PEZCG,pi);

    TACCM_CGAA_MD_Propagetor_NH_MP1998_Z_2(Z,velZ,massZ,thetaAA,thetaCG,zetaZ,V_zetaZ,
					   QZ,NfKTZ,numZ,&KEZ,&KEvZ,&PEvZ,
					   dt,dt2,nc,wdt4,wdt2,KZAA,KZCG,
					   PEZAA,PEZCG,PEZ,frcZ,pi);
      
    if (i%interval==0) {
      KEAA=KEAA/UNITT;     TAA=KEAA/((3*numatom)*k_B)*2.0;
      PEvAA=PEvAA/UNITT;   KEvAA=KEvAA/UNITT;

      KECG=KECG/UNITT;     TCG=KECG/((3*/*numatom*/numheavyatom)*k_B)*2.0;
      PEvCG=PEvCG/UNITT;  KEvCG=KEvCG/UNITT;

      PEAA=0.5*e.p_e_t+0.5*e.p_LJ_t+0.5*e.p_e_14_t+0.5*e.p_LJ_14_t+e.p_d_t+e.p_a_t+e.p_b_t;
      EtAA=PEAA+KEAA+PEvAA+KEvAA;
      fprintf(outputfileAA,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "
                  	   ,i+1,PEAA, KEAA, KEvAA,PEvAA,EtAA,*PEZAA,TAA);

      PECG=e_CG.p_t;
      EtCG=PECG+KECG+PEvCG+KEvCG;
      fprintf(outputfileCG,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "
  	                   ,i+1,PECG,KECG, KEvCG, PEvCG,EtCG,*PEZCG,TCG);

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
      fprintf(outputfileAA,"%d %e %e %e %e %e %e \n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);

      fprintf(outputfileCG,"%d %e %e %e %e %e %e \n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);

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


