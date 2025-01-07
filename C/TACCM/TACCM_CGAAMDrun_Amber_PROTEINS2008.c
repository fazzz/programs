
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "FFL.h"

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

#include "TACCM_CGAAMDrun_Amber_PROTEINS2008.h"
//#include "TACCM_CGAAMDrun_test_CG.h"
#include "TACCM_MDrun.h"
#include "TACCM_MD.h"

#include "MD_NHC_MP1996.h"
#include "MDrun.h"
#include "MD.h"

double runTACCM_CGAA_MD_NHC_MP1998_Amber_PROTEINS2008(// AA ////////////////////////////////////////////
						      double *crdAA,double *velAA, 
						      double *zetaAA,double *V_zetaAA, double QAA,
						      struct potential e, struct force f, double TAA, double NfKTAA,
						      double *avePEAA, double *aveKEAA,double *aveTAA,
						      double *varPEAA, double *varKEAA,double *varTAA,
						      struct my_netcdf_out_id_MCD nc_id_MCDAA,  FILE *outputfileAA,
						      // CG ////////////////////////////////////////////
						      double *crdCG,double *velCG, 
						      double *zetaCG,double *V_zetaCG, double QCG,
						      struct potential_GOLMAA_PROTEINS2008 e_CG, 
						      double TCG, double NfKTCG,
						      double *avePECG, double *aveKECG,double *aveTCG,
						      double *varPECG, double *varKECG,double *varTCG,
						      struct my_netcdf_out_id_MCD nc_id_MCDCG,  FILE *outputfileCG,
						      // Z  /////////////////////////////////////////////
						      double *Z,double *velZ,double massZ,
						      double *zetaZ,double *V_zetaZ,
						      double QZ,double NfKTZ,double TZ,
						      int numZ,double KZAA, double KZCG,int **pairs,
						      double *avePEZ, double *aveKEZ,double *aveTZ,
						      double *varPEZ, double *varKEZ,double *varTZ, 
						      FILE *trjfileZ, FILE *trjfilThetaAA, FILE *trjfilThetaCG,
						      // CM  ///////////////////////////////////////////////
						      double *mass, int numatom, int numheavyatom,
						      int numstep, int interval,int *l,
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

  ffL_calcffandforce(crdAA,numatom,&e,&f);
  //  ffL_calcffandforce_CG(crdCG,numatom,&e_CG,&f_CG,parameterCG);
  GOLMAA_PROTEINS2008_ff_calcff_b(crdCG,numatom,&e_CG);

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

    PECG=TACCM_MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008(crdCG,velCG,mass,zetaCG,V_zetaCG,
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

double runTACCM_CGAA_MD_NHC_MP1998_Amber_PROTEINS2008_Z_asCA(// AA ////////////////////////////////////////////
							     double *crdAA,double *velAA, 
							     double *zetaAA,double *V_zetaAA, double QAA,
							     struct potential e, struct force f, double TAA, double NfKTAA,
							     double *avePEAA, double *aveKEAA,double *aveTAA,
							     double *varPEAA, double *varKEAA,double *varTAA,
							     struct my_netcdf_out_id_MCD nc_id_MCDAA,  FILE *outputfileAA,
							     // CG ////////////////////////////////////////////
							     double *crdCG,double *velCG, 
							     double *zetaCG,double *V_zetaCG, double QCG,
							     struct potential_GOLMAA_PROTEINS2008 e_CG, 
							     double TCG, double NfKTCG,
							     double *avePECG, double *aveKECG,double *aveTCG,
							     double *varPECG, double *varKECG,double *varTCG,
							     struct my_netcdf_out_id_MCD nc_id_MCDCG,  FILE *outputfileCG,
							     // Z  /////////////////////////////////////////////
							     double *Z,double *velZ,double massZ,
							     double *zetaZ,double *V_zetaZ,
							     double QZ,double NfKTZ,double TZ,
							     int numZ,double KZAA, double KZCG,int *index,
							     double *avePEZ, double *aveKEZ,double *aveTZ,
							     double *varPEZ, double *varKEZ,double *varTZ, 
							     struct my_netcdf_out_id_MCD trjfileZ, 
							     struct my_netcdf_out_id_MCD trjfilThetaAA,
							     struct my_netcdf_out_id_MCD trjfilThetaCG,
							     // CM  ///////////////////////////////////////////////
							     double *mass, int numatom, int numheavyatom,
							     int numstep, int interval,int *l,
							     double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
							     double UNITT, double k_B,double pi,
							     double *PEZAA, double *PEZCG, double *PEZ) {
  int i,j,k,ll;

  double PEAA=0.0,KEAA=0.0,EtAA,PEvAA,KEvAA;
  double PECG=0.0,KECG=0.0,EtCG,PEvCG,KEvCG;
  double KEZ,KEvZ,PEvZ,EtZ;

  double *thetaAA,*thetaCG,*frcZ,*fAA,*fCG;
  double summass,COM[3],crd_nc[MAXATOM][3],crd_ncZ[MAXATOM][3];

  ffL_calcffandforce(crdAA,numatom,&e,&f);
  //  ffL_calcffandforce_CG(crdCG,numatom,&e_CG,&f_CG,parameterCG);
  GOLMAA_PROTEINS2008_ff_calcff_b(crdCG,numatom,&e_CG);

  fAA=(double *)gcemalloc(sizeof(double)*numZ*3);
  fCG=(double *)gcemalloc(sizeof(double)*numZ*3);
  frcZ=(double *)gcemalloc(sizeof(double)*numZ*3);

  thetaAA=(double *)gcemalloc(sizeof(double)*numZ*3);
  thetaCG=(double *)gcemalloc(sizeof(double)*numZ*3);

  k=0;
  for (i=0;i<numatom;++i) {
    if (i+1==index[k]) {
      for (j=0;j<3;++j) {
	thetaAA[k*3+j]=crdAA[i*3+j];
	thetaCG[k*3+j]=crdCG[i*3+j];
      }
      ++k;
    }
  }

  for (i=0;i<numZ;++i) {
    for (j=0;j<3;++j) {
      fAA[/*(index[*/i/*]-1)*/*3+j]=-KZAA*(thetaAA[i*3+j]-Z[i*3+j])*418.4070;
      fCG[/*(index[*/i/*]-1)*/*3+j]=-KZCG*(thetaCG[i*3+j]-Z[i*3+j])*418.4070;
    }
  }

  for (i=0;i<numZ*3;++i) frcZ[i]=fAA[i]+fCG[i];

  for (i=0;i<numstep;++i) {
    k=0;
    for (j=0;j<numatom;++j) {
      if (j+1==index[k]) {
	for (ll=0;ll<3;++ll) {
	  thetaAA[k*3+ll]=crdAA[j*3+ll];
	  thetaCG[k*3+ll]=crdCG[j*3+ll];
	}
	++k;
      }
    }

    PEAA=TACCM_MD_Propagetor_NH_MP1998_AAFF_Amber_Z_asCA(crdAA,velAA,mass,zetaAA,V_zetaAA,
							 QAA,NfKTAA,numatom,&KEAA,&KEvAA,&PEvAA,
							 dt,dt2,nc,wdt4,wdt2,
							 &e,&f,Z,numZ,thetaAA,KZAA,index,PEZAA,pi);

    PECG=TACCM_MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008_Z_asCA(crdCG,velCG,mass,zetaCG,V_zetaCG,
    								  QCG,NfKTCG,numatom,&KECG,&KEvCG,&PEvCG,
    								  dt,dt2,nc,wdt4,wdt2,
    								  &e_CG,Z,numZ,thetaCG,KZCG,index,PEZCG,pi);

    TACCM_CGAA_MD_Propagetor_NH_MP1998_Z_2_asCA(Z,velZ,massZ,thetaAA,thetaCG,zetaZ,V_zetaZ,
    						QZ,NfKTZ,numZ,&KEZ,&KEvZ,&PEvZ,
    						dt,dt2,nc,wdt4,wdt2,
    						KZAA,KZCG,PEZAA,PEZCG,PEZ,frcZ);
      
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
      //      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) crdAA[j*3+k]-=COM[k];
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crdAA[j*3+k];
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) crd_nc[j][k]-=COM[k];
      myncL_put_crd_ene_MCD(nc_id_MCDAA,*l,crd_nc,e,0.0);

      myncL_put_crd_ene_MCD(nc_id_MCDCG,*l,crd_nc,e,0.0);

      summass=0.0; for (j=0;j<numZ;++j) summass+=massZ;
      for (j=0;j<3;++j) COM[j]=0.0;
      for (j=0;j<numZ;++j)  for (k=0;k<3;++k) COM[k]+=Z[j*3+k]*massZ/summass;
      /**************************************/
      /* for (j=0;j<numZ;++j) {		    */
      /* 	for (k=0;k<3;++k) {	    */
      /* 	  Z[j*3+k]-=COM[k];	    */
      /* 	  thetaAA[j*3+k]-=COM[k];   */
      /* 	  thetaCG[j*3+k]-=COM[k];   */
      /* 	}			    */
      /* }				    */
      /**************************************/
      for (j=0;j<numZ;++j) for (k=0;k<3;++k) crd_ncZ[j][k]=Z[j*3+k]-COM[k];
      myncL_put_crd_ene_MCD(trjfileZ,*l,crd_ncZ,e,0.0);
      for (j=0;j<numZ;++j) for (k=0;k<3;++k) crd_ncZ[j][k]=thetaAA[j*3+k]-COM[k];
      myncL_put_crd_ene_MCD(trjfilThetaAA,*l,crd_ncZ,e,0.0);
      for (j=0;j<numZ;++j) for (k=0;k<3;++k) crd_ncZ[j][k]=thetaCG[j*3+k]-COM[k];
      myncL_put_crd_ene_MCD(trjfilThetaCG,*l,crd_ncZ,e,0.0);

      ++(*l);

      ///////////////// TACCM //////////////////////
      KEZ=KEZ/UNITT;      TZ=KEZ/(numZ*3*k_B)*2.0;

      PEvZ=PEvZ/UNITT;      KEvZ=KEvZ/UNITT;

      EtZ=*PEZ+KEZ+PEvZ+KEvZ;
      fprintf(outputfileAA,"%d %e %e %e %e %e %e \n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);

      fprintf(outputfileCG,"%d %e %e %e %e %e %e \n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);

      ///////////////// TACCM //////////////////////
    }
  }

  return *PEZ;
}

double runTACCM_CGAA_MD_NHC_MP1998_2PROTEINS2008(// CG1 ////////////////////////////////////////////
						 double *crdCG1,double *velCG1, 
						 double *zetaCG1,double *V_zetaCG1, double QCG1,
						 struct potential_GOLMAA_PROTEINS2008 e_CG1, 
						 double TCG1, double NfKTCG1,
						 double *avePECG1, double *aveKECG1,double *aveTCG1,
						 double *varPECG1, double *varKECG1,double *varTCG1,
						 struct my_netcdf_out_id_MCD nc_id_MCDCG1,  FILE *outputfileCG1,
						 // CG2 ////////////////////////////////////////////
						 double *crdCG2,double *velCG2, 
						 double *zetaCG2,double *V_zetaCG2, double QCG2,
						 struct potential_GOLMAA_PROTEINS2008 e_CG2, 
						 double TCG2, double NfKTCG2,
						 double *avePECG2, double *aveKECG2,double *aveTCG2,
						 double *varPECG2, double *varKECG2,double *varTCG2,
						 struct my_netcdf_out_id_MCD nc_id_MCDCG2,  FILE *outputfileCG2,
						 // Z  /////////////////////////////////////////////
						 double *Z,double *velZ,double massZ,
						 double *zetaZ,double *V_zetaZ,
						 double QZ,double NfKTZ,double TZ,
						 int numZ,double KZCG1, double KZCG2,int **pairs,
						 double *avePEZ, double *aveKEZ,double *aveTZ,
						 double *varPEZ, double *varKEZ,double *varTZ, 
						 FILE *trjfileZ, FILE *trjfilThetaCG1, FILE *trjfilThetaCG2,
						 // CM  ///////////////////////////////////////////////
						 double *mass, int numatom, int numheavyatom,
						 int numstep, int interval,int *l,
						 double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
						 double UNITT, double k_B,double pi,
						 double *PEZCG1, double *PEZCG2, double *PEZ) {
  int i,j,k;

  double PECG1=0.0,KECG1=0.0,EtCG1,PEvCG1,KEvCG1;
  double PECG2=0.0,KECG2=0.0,EtCG2,PEvCG2,KEvCG2;
  double KEZ,KEvZ,PEvZ,EtZ;

  double *thetaCG1,*thetaCG2,*frcZ,*fCG1,*fCG2;
  double summass,COM[3],crd_nc[MAXATOM][3];

  struct potential e;

  GOLMAA_PROTEINS2008_ff_calcff_b(crdCG1,numatom,&e_CG1);
  GOLMAA_PROTEINS2008_ff_calcff_b(crdCG2,numatom,&e_CG2);

  fCG1=(double *)gcemalloc(sizeof(double)*numZ);
  fCG2=(double *)gcemalloc(sizeof(double)*numZ);
  frcZ=(double *)gcemalloc(sizeof(double)*numZ);

  thetaCG1=(double *)gcemalloc(sizeof(double)*numZ);
  thetaCG2=(double *)gcemalloc(sizeof(double)*numZ);

  TACCM_CTheta(crdCG1,numatom,thetaCG1,numZ,pairs,pi);
  TACCM_CTheta(crdCG2,numatom,thetaCG2,numZ,pairs,pi);
  TACCM_calc_eff_FF_Z(Z,numZ,thetaCG1,KZCG1,fCG1,pi);
  TACCM_calc_eff_FF_Z(Z,numZ,thetaCG2,KZCG2,fCG2,pi);

  for (i=0;i<numZ;++i) frcZ[i]=fCG1[i]+fCG2[i];

  for (i=0;i<numstep;++i) {
    TACCM_CTheta(crdCG1,numatom,thetaCG1,numZ,pairs,pi);
    TACCM_CTheta(crdCG2,numatom,thetaCG2,numZ,pairs,pi);

    PECG1=TACCM_MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008(crdCG1,velCG1,mass,zetaCG1,V_zetaCG1,
							   QCG1,NfKTCG1,numatom,&KECG1,&KEvCG1,&PEvCG1,
							   dt,dt2,nc,wdt4,wdt2,
							   &e_CG1,Z,numZ,thetaCG1,KZCG1,pairs,PEZCG1,pi);


    PECG2=TACCM_MD_Propagetor_NH_MP1998_GOLMAA_PROTEINS2008(crdCG2,velCG2,mass,zetaCG2,V_zetaCG2,
							   QCG2,NfKTCG2,numatom,&KECG2,&KEvCG2,&PEvCG2,
							   dt,dt2,nc,wdt4,wdt2,
							   &e_CG2,Z,numZ,thetaCG2,KZCG2,pairs,PEZCG2,pi);

    TACCM_CGAA_MD_Propagetor_NH_MP1998_Z_2(Z,velZ,massZ,thetaCG1,thetaCG2,zetaZ,V_zetaZ,
					   QZ,NfKTZ,numZ,&KEZ,&KEvZ,&PEvZ,
					   dt,dt2,nc,wdt4,wdt2,KZCG1,KZCG2,
					   PEZCG1,PEZCG2,PEZ,frcZ,pi);
      
    if (i%interval==0) {
      KECG1=KECG1/UNITT;     TCG1=KECG1/((3*/*numatom*/numheavyatom)*k_B)*2.0;
      PEvCG1=PEvCG1/UNITT;  KEvCG1=KEvCG1/UNITT;

      KECG2=KECG2/UNITT;     TCG2=KECG2/((3*/*numatom*/numheavyatom)*k_B)*2.0;
      PEvCG2=PEvCG2/UNITT;  KEvCG2=KEvCG2/UNITT;

      PECG1=e_CG1.p_t;
      EtCG1=PECG1+KECG1+PEvCG1+KEvCG1;
      fprintf(outputfileCG1,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "
                  	   ,i+1,PECG1, KECG1, KEvCG1,PEvCG1,EtCG1,*PEZCG1,TCG1);

      PECG2=e_CG2.p_t;
      EtCG2=PECG2+KECG2+PEvCG2+KEvCG2;
      fprintf(outputfileCG2,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "
  	                   ,i+1,PECG2,KECG2, KEvCG2, PEvCG2,EtCG2,*PEZCG2,TCG2);

      //      *avePECG1=(i*(*avePECG1)+PECG1)/(i+1); *varPECG1=(i*(*varPECG1)+PECG1*PECG1)/(i+1);
      //      *aveKECG1=(i*(*aveKECG1)+KECG1)/(i+1); *varKECG1=(i*(*varKECG1)+KECG1*KECG1)/(i+1);
      //      *aveTCG1=(i*(*aveTCG1)+TCG1)/(i+1);  *varTCG1=(i*(*varTCG1)+TCG1*TCG1)/(i+1);

      //      *avePECG2=(i*(*avePECG2)+PECG2)/(i+1); *varPECG2=(i*(*varPECG2)+PECG2*PECG2)/(i+1);
      //      *aveKECG2=(i*(*aveKECG2)+KECG2)/(i+1); *varKECG2=(i*(*varKECG2)+KECG2*KECG2)/(i+1);
      //      *aveTCG2=(i*(*aveTCG2)+TCG2)/(i+1);  *varTCG2=(i*(*varTCG2)+TCG2*TCG2)/(i+1);

      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crdCG1[j*3+k];
      myncL_put_crd_ene_MCD(nc_id_MCDCG1,*l,crd_nc,e,0.0);

      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crdCG2[j*3+k];
      myncL_put_crd_ene_MCD(nc_id_MCDCG2,*l,crd_nc,e,0.0);
      ++(*l);

      ///////////////// TACCM //////////////////////
      KEZ=KEZ/UNITT;      TZ=KEZ/(numZ*k_B)*2.0;

      PEvZ=PEvZ/UNITT;      KEvZ=KEvZ/UNITT;

      EtZ=*PEZ+KEZ+PEvZ+KEvZ;
      fprintf(outputfileCG1,"%d %e %e %e %e %e %e \n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);

      fprintf(outputfileCG2,"%d %e %e %e %e %e %e \n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);

      for (j=0;j<numZ;++j) fprintf(trjfileZ,"%e ",Z[j]);
      fprintf(trjfileZ,"\n");
      for (j=0;j<numZ;++j) fprintf(trjfilThetaCG1,"%e ",thetaCG1[j]);
      fprintf(trjfilThetaCG1,"\n");      
      for (j=0;j<numZ;++j) fprintf(trjfilThetaCG2,"%e ",thetaCG2[j]);
      fprintf(trjfilThetaCG2,"\n");      
      ///////////////// TACCM //////////////////////
    }
  }

  return *PEZ;
}
