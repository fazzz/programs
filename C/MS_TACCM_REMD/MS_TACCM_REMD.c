
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"

#include "MD_NHC_MP1996.h"
#include "MDrun.h"
#include "MD.h"

double run_multi_MS_TACCM_Amber_PROTEINS2008_pca(int numAA, int numCG,
						 struct AADataMSREMD_Amber *FGdata,
						 struct CGDataMSREMD_PROTEINS2008 *CGdata,
						 struct ZDataMSREMD_A_P2008_pca Zdata,
						 struct AmberParmL *ap,struct AmberParmL *ap_CG,
						 double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
						 double UNITT, double k_B,double pi,
						 double *PEZAA, double *PEZCG,double *PEZ){
  int i,j,k;

  double *PEAA,*KEAA,*EtAA,*PEvAA,*KEvAA;
  double *PECG,*KECG,*EtCG,*PEvCG,*KEvCG;

  double KEZ,KEvZ,PEvZ,EtZ;

  double **thetaAA,**thetaCG,**thetaCG;
  double *frcZ,**fAA,**fCG;

  double summass,COM[3],crd_nc[MAXATOM][3];

  fAA=(double **)gcemalloc(sizeof(double *)*numAA);
  thetaAA=(double **)gcemalloc(sizeof(double *)*numAA);
  for (i=0;i<numAA;++i) {
    ffL_calcffandforce(FGdata[i].crd,ap[i].NATOM,&(FGdata[i].e),&f(FGdata[i].f));
    fAA[i]=(double *)gcemalloc(sizeof(double)*numZ);
    thetaAA[i]=(double *)gcemalloc(sizeof(double)*numZ);
    TACCM_pca_CTheta(FGdata[i].crd,ap[i].NATOM,thetaAA[i],numZ,pairs,pi);
  }

  fCG=(double **)gcemalloc(sizeof(double *)*numCG);
  thetaCCG=(double **)gcemalloc(sizeof(double *)*numCG);
  for (i=0;i<numCG;++i) {
    GOLMAA_PROTEINS2008_ff_calcff_b(CGdata[i].crd,ap_CG[i].NATOM,&(CGdata[i].e_CG));
    fCG[i]=(double *)gcemalloc(sizeof(double)*numZ);
    thetaCG[i]=(double *)gcemalloc(sizeof(double)*numZ);
  }

  frcZ=(double *)gcemalloc(sizeof(double)*numZ);

  TACCM_pca_CTheta(crdCG1,numatom,thetaCG1,numZ,pairs,pi);

  TACCM_calc_eff_FF_Z(Z,numZ,thetaAA,KZAA,fAA,pi);
  TACCM_calc_eff_FF_Z(Z,numZ,thetaCG1,KZCG1,fCG1,pi);
  TACCM_calc_eff_FF_Z(Z,numZ,thetaCG2,KZCG2,fCG2,pi);

  for (i=0;i<numZ;++i) {
    frcZ[i]=fAA[i]+fCG1[i]+fCG2[i];
  }

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
