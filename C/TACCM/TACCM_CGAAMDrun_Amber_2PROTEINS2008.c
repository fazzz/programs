
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "FFL.h"

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

#include "TACCM_CGAAMDrun_Amber_PROTEINS2008.h"
#include "TACCM_CGAAMDrun_Amber_2PROTEINS2008.h"
//#include "TACCM_CGAAMDrun_test_CG.h"
#include "TACCM_MDrun.h"
#include "TACCM_MD.h"

#include "MD_NHC_MP1996.h"
#include "MDrun.h"
#include "MD.h"

double runTACCM_2CG1FG_MD_NHC_MP1998_Amber_PROTEINS2008(// AA ////////////////////////////////////////////
							double *crdAA,double *velAA, 
							double *zetaAA,double *V_zetaAA, double QAA,
							struct potential e, struct force f, double TAA, double NfKTAA,
							double *avePEAA, double *aveKEAA,double *aveTAA,
							double *varPEAA, double *varKEAA,double *varTAA,
							struct my_netcdf_out_id_MCD nc_id_MCDAA,  FILE *outputfileAA,
							// CG1 ///////////////////////////////////////////
							double *crdCG1,double *velCG1, 
							double *zetaCG1,double *V_zetaCG1, double QCG1,
							struct potential_GOLMAA_PROTEINS2008 e_CG1, 
							double TCG1, double NfKTCG1,
							double *avePECG1, double *aveKECG1,double *aveTCG1,
							double *varPECG1, double *varKECG1,double *varTCG1,
							struct my_netcdf_out_id_MCD nc_id_MCDCG1,
							FILE *outputfileCG1,
							// CG2 ///////////////////////////////////////////
							double *crdCG2,double *velCG2, 
							double *zetaCG2,double *V_zetaCG2, double QCG2,
							struct potential_GOLMAA_PROTEINS2008 e_CG2, 
							double TCG2, double NfKTCG2,
							double *avePECG2, double *aveKECG2,double *aveTCG2,
							double *varPECG2, double *varKECG2,double *varTCG2,
							struct my_netcdf_out_id_MCD nc_id_MCDCG2,
							FILE *outputfileCG2,
							// Z  /////////////////////////////////////////////
							double *Z,double *velZ,double massZ,
							double *zetaZ,double *V_zetaZ,
							double QZ,double NfKTZ,double TZ,
							int numZ,double KZAA, double KZCG1,double KZCG2,int **pairs,
							double *avePEZ, double *aveKEZ,double *aveTZ,
							double *varPEZ, double *varKEZ,double *varTZ, 
							FILE *trjfileZ, FILE *trjfilThetaAA,
							FILE *trjfilThetaCG1,FILE *trjfilThetaCG2,
							// CM  ///////////////////////////////////////////////
							double *mass, int numatom, int numheavyatom,
							int numstep, int interval,int *l,
							double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
							double UNITT, double k_B,double pi,
							double *PEZAA, double *PEZCG1, double *PEZCG2,
							double *PEZ) {
  int i,j,k;

  double PEAA=0.0,KEAA=0.0,EtAA,PEvAA,KEvAA;
  double PECG1=0.0,KECG1=0.0,EtCG1,PEvCG1,KEvCG1;
  double PECG2=0.0,KECG2=0.0,EtCG2,PEvCG2,KEvCG2;
  double KEZ,KEvZ,PEvZ,EtZ;

  double *thetaAA,*thetaCG1,*thetaCG2,*frcZ,*fAA,*fCG1,*fCG2;
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
  GOLMAA_PROTEINS2008_ff_calcff_b(crdCG1,numatom,&e_CG1);
  GOLMAA_PROTEINS2008_ff_calcff_b(crdCG2,numatom,&e_CG2);

  fAA=(double *)gcemalloc(sizeof(double)*numZ);
  fCG1=(double *)gcemalloc(sizeof(double)*numZ);
  fCG2=(double *)gcemalloc(sizeof(double)*numZ);
  frcZ=(double *)gcemalloc(sizeof(double)*numZ);

  thetaAA=(double *)gcemalloc(sizeof(double)*numZ);
  thetaCG1=(double *)gcemalloc(sizeof(double)*numZ);
  thetaCG2=(double *)gcemalloc(sizeof(double)*numZ);

  TACCM_CTheta(crdAA,numatom,thetaAA,numZ,pairs,pi);
  TACCM_CTheta(crdCG1,numatom,thetaCG1,numZ,pairs,pi);
  TACCM_CTheta(crdCG2,numatom,thetaCG2,numZ,pairs,pi);
  TACCM_calc_eff_FF_Z(Z,numZ,thetaAA,KZAA,fAA,pi);
  TACCM_calc_eff_FF_Z(Z,numZ,thetaCG1,KZCG1,fCG1,pi);
  TACCM_calc_eff_FF_Z(Z,numZ,thetaCG2,KZCG2,fCG2,pi);

  for (i=0;i<numZ;++i) frcZ[i]=fAA[i]+fCG1[i]+fCG2[i];

  for (i=0;i<numstep;++i) {
    TACCM_CTheta(crdAA,numatom,thetaAA,numZ,pairs,pi);
    TACCM_CTheta(crdCG1,numatom,thetaCG1,numZ,pairs,pi);
    TACCM_CTheta(crdCG2,numatom,thetaCG2,numZ,pairs,pi);

    PEAA=TACCM_MD_Propagetor_NH_MP1998_AAFF_Amber(crdAA,velAA,mass,zetaAA,V_zetaAA,
						  QAA,NfKTAA,numatom,&KEAA,&KEvAA,&PEvAA,
						  dt,dt2,nc,wdt4,wdt2,
						  &e,&f,Z,numZ,thetaAA,KZAA,pairs,PEZAA,pi);

    PECG1=TACCM_MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008(crdCG1,velCG1,mass,zetaCG1,V_zetaCG1,
							    QCG1,NfKTCG1,numatom,&KECG1,&KEvCG1,&PEvCG1,
							    dt,dt2,nc,wdt4,wdt2,
							    &e_CG1,Z,numZ,thetaCG1,KZCG1,pairs,PEZCG1,pi);

    PECG2=TACCM_MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008(crdCG2,velCG2,mass,zetaCG2,V_zetaCG2,
							    QCG2,NfKTCG2,numatom,&KECG2,&KEvCG2,&PEvCG2,
							    dt,dt2,nc,wdt4,wdt2,
							    &e_CG2,Z,numZ,thetaCG2,KZCG2,pairs,PEZCG2,pi);

    TACCM_2CG1FG_MD_Propagetor_NH_MP1998_Z_2(Z,velZ,massZ,thetaAA,thetaCG1,thetaCG2,zetaZ,V_zetaZ,
					     QZ,NfKTZ,numZ,&KEZ,&KEvZ,&PEvZ,
					     dt,dt2,nc,wdt4,wdt2,KZAA,KZCG1,KZCG2,
					     PEZAA,PEZCG1,PEZCG2,PEZ,frcZ,pi);
      
    if (i%interval==0) {
      KEAA=KEAA/UNITT;     TAA=KEAA/((3*numatom)*k_B)*2.0;
      PEvAA=PEvAA/UNITT;   KEvAA=KEvAA/UNITT;

      KECG1=KECG1/UNITT;     TCG1=KECG1/((3*/*numatom*/numheavyatom)*k_B)*2.0;
      PEvCG1=PEvCG1/UNITT;  KEvCG1=KEvCG1/UNITT;

      KECG2=KECG2/UNITT;     TCG2=KECG2/((3*/*numatom*/numheavyatom)*k_B)*2.0;
      PEvCG2=PEvCG2/UNITT;  KEvCG2=KEvCG2/UNITT;

      PEAA=0.5*e.p_e_t+0.5*e.p_LJ_t+0.5*e.p_e_14_t+0.5*e.p_LJ_14_t+e.p_d_t+e.p_a_t+e.p_b_t;
      EtAA=PEAA+KEAA+PEvAA+KEvAA;
      fprintf(outputfileAA,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "
                  	   ,i+1,PEAA, KEAA, KEvAA,PEvAA,EtAA,*PEZAA,TAA);

      PECG1=e_CG1.p_t;
      EtCG1=PECG1+KECG1+PEvCG1+KEvCG1;
      fprintf(outputfileCG1,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "
  	                   ,i+1,PECG1,KECG1, KEvCG1, PEvCG1,EtCG1,*PEZCG1,TCG1);

      PECG2=e_CG2.p_t;
      EtCG2=PECG2+KECG2+PEvCG2+KEvCG2;
      fprintf(outputfileCG2,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "
  	                   ,i+1,PECG2,KECG2, KEvCG2, PEvCG2,EtCG2,*PEZCG2,TCG2);

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

      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crdCG1[j*3+k];
      myncL_put_crd_ene_MCD(nc_id_MCDCG1,*l,crd_nc,e,0.0);

      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crdCG2[j*3+k];
      myncL_put_crd_ene_MCD(nc_id_MCDCG2,*l,crd_nc,e,0.0);
      ++(*l);

      ///////////////// TACCM //////////////////////
      KEZ=KEZ/UNITT;      TZ=KEZ/(numZ*k_B)*2.0;

      PEvZ=PEvZ/UNITT;      KEvZ=KEvZ/UNITT;

      EtZ=*PEZ+KEZ+PEvZ+KEvZ;
      fprintf(outputfileAA,"%d %e %e %e %e %e %e \n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);

      fprintf(outputfileCG1,"%d %e %e %e %e %e %e \n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);

      fprintf(outputfileCG2,"%d %e %e %e %e %e %e \n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);

      for (j=0;j<numZ;++j) fprintf(trjfileZ,"%e ",Z[j]);
      fprintf(trjfileZ,"\n");
      for (j=0;j<numZ;++j) fprintf(trjfilThetaAA,"%e ",thetaAA[j]);
      fprintf(trjfilThetaAA,"\n");      
      for (j=0;j<numZ;++j) fprintf(trjfilThetaCG1,"%e ",thetaCG1[j]);
      fprintf(trjfilThetaCG1,"\n");      
      for (j=0;j<numZ;++j) fprintf(trjfilThetaCG2,"%e ",thetaCG2[j]);
      fprintf(trjfilThetaCG2,"\n");      
      ///////////////// TACCM //////////////////////
    }
  }

  return *PEZ;
}

double TACCM_2CG1FG_MD_Propagetor_NH_MP1998_Z_2(double *Z,double *velZ,double massZ,
						double *thetaAA,double *thetaCG1,double *thetaCG2,
						double *zeta,double *V_zeta,double Q,double NfKT,int numZ,
						double *KE, double *KEv,double *PEv,
						double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
						double KZAA, double KZCG1, double KZCG2,
						double *PEZAA,double *PEZCG1,double *PEZCG2,
						double *PEZ,double *f,double pi) {
  int i,j,k;
  //  double PEZAA,PEZCG;
  double *fAA,*fCG1,*fCG2;

  fAA=(double *)gcemalloc(sizeof(double)*numZ);
  fCG1=(double *)gcemalloc(sizeof(double)*numZ);
  fCG2=(double *)gcemalloc(sizeof(double)*numZ);

  TACCM_MD_Propagetor_NH_Single_part_MP1996(velZ,massZ,zeta,V_zeta,Q,NfKT,numZ,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numZ;++i) velZ[i]+=dt2*f[i]/massZ;
  for (i=0;i<numZ;++i) {
    Z[i]+=dt*velZ[i];
    while (Z[i]>pi) Z[i]-=2.0*pi;
    while (Z[i]<=-pi) Z[i]+=2.0*pi;
  }
  TACCM_calc_eff_FF_Z_2(Z,numZ,thetaAA,KZAA,fAA,PEZAA,pi);
  TACCM_calc_eff_FF_Z_2(Z,numZ,thetaCG1,KZCG1,fCG1,PEZCG1,pi);
  TACCM_calc_eff_FF_Z_2(Z,numZ,thetaCG2,KZCG2,fCG2,PEZCG2,pi);

  *PEZ=(*PEZAA)+(*PEZCG1)+(*PEZCG2); for (i=0;i<numZ;++i) f[i]=fAA[i]+fCG1[i]+fCG2[i];

  for (i=0;i<numZ;++i) velZ[i]+=dt2*f[i]/massZ;
  //////////////////////////////////////////////////////////////////////////////////////

  TACCM_MD_Propagetor_NH_Single_part_MP1996(velZ,massZ,zeta,V_zeta,Q,NfKT,numZ,nc,wdt4,wdt2);

  *KE=0.0; for (i=0;i<numZ;++i) *KE+=0.5*massZ*velZ[i]*velZ[i];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return 0.0;
}
