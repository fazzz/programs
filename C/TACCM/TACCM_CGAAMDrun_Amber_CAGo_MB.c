
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <math.h>

#include "EF.h"
#include "FFL.h"
#include "PTL.h"

#include "TOPO.h"
#include "LA.h"

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

#include "GOLM_Clementi_MB_wmCutOff.h"

#include "TACCM_CGAAMDrun_Amber_CAGo_MB.h"
#include "TACCM_CGAAMDrun_Amber_CAGo_MB_CTheta.h"

#include "TACCM_CGAAMDrun_test.h"
#include "TACCM_CGAAMDrun_test_CG.h"
#include "TACCM_MDrun.h"
#include "TACCM_MD.h"

#include "MD_NHC_MP1996.h"
#include "MDrun.h"
#include "MD.h"

#include "FVDIHED.h"

//#include "MuSTARMD_hs_funcs.h"

double TACCM_CGAA_MD_Propagetor_NH_MP1998_Z_Amber_CAGo_MB_2(double *Z,double *velZ,double massZ,
							    double *thetaAA,double *thetaCG,
							    double *zeta,double *V_zeta,double Q,double NfKT,
							    int numZ_dih, int numZ_ang, int numZ_bon,
							    double *KE, double *KEv,double *PEv,
							    double dt,double dt2,int nc,
							    double wdt4[3],double wdt2[3],
							    double KZAA, double KZCG,
							    double *PEZAA,double *PEZCG,double *PEZ,double *f,
							    double UNITT, double pi) {
  int i,j,k;
  //  double PEZAA,PEZCG;
  double *fAA,*fCG;

  fAA=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));
  fCG=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));

  TACCM_MD_Propagetor_NH_Single_part_MP1996(velZ,massZ,zeta,V_zeta,Q,NfKT,
					    numZ_dih+numZ_ang+numZ_bon,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numZ_dih+numZ_ang+numZ_bon;++i) velZ[i]+=dt2*f[i]/massZ;
  for (i=0;i<numZ_dih+numZ_ang+numZ_bon;++i) {
    Z[i]+=dt*velZ[i];
    while (Z[i]>pi) Z[i]-=2.0*pi;
    while (Z[i]<=-pi) Z[i]+=2.0*pi;
  }
  TACCM_calc_eff_FF_Z_2_Amber_CAGo_MB(Z,numZ_dih,numZ_ang,numZ_bon,
				      thetaAA,KZAA,fAA,PEZAA,UNITT,pi);
  TACCM_calc_eff_FF_Z_2_Amber_CAGo_MB(Z,numZ_dih,numZ_ang,numZ_bon,
				      thetaCG,KZCG,fCG,PEZCG,UNITT,pi);

  *PEZ=(*PEZAA)+(*PEZCG); for (i=0;i<numZ_dih+numZ_ang+numZ_bon;++i) f[i]=fAA[i]+fCG[i];

  for (i=0;i<numZ_dih+numZ_ang+numZ_bon;++i) velZ[i]+=dt2*f[i]/massZ;
  //////////////////////////////////////////////////////////////////////////////////////

  TACCM_MD_Propagetor_NH_Single_part_MP1996(velZ,massZ,zeta,V_zeta,Q,NfKT,
					    numZ_dih+numZ_ang+numZ_bon,nc,wdt4,wdt2);

  *KE=0.0; for (i=0;i<numZ_dih+numZ_ang+numZ_bon;++i) *KE+=0.5*massZ*velZ[i]*velZ[i];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return 0.0;
}

double runTACCM_CGAA_MD_NHC_MP1998_Amber_CAGo_MB_2(// AA /////////////////////////////////////////////////////////
						   double *crdAA,double *velAA, 
						   double *zetaAA,double *V_zetaAA, double QAA,
						   struct potential e, struct force f, struct AmberParmL ap_AA,
						   double TAA, double NfKTAA,
						   double *avePEAA, double *aveKEAA,double *aveTAA,
						   double *varPEAA, double *varKEAA,double *varTAA,
						   struct my_netcdf_out_id_MCD nc_id_MCDAA,  FILE *outputfileAA,
						   // CG /////////////////////////////////////////////////////////
						   double *crdCG,double *velCG, 
						   double *zetaCG,double *V_zetaCG, double QCG,
						   struct potential_GOLM_Clementi_MB e_CG,
						   double de, double d2,
						   struct AmberParmL ap_CG,
						   double TCG, double NfKTCG,
						   double *avePECG, double *aveKECG,double *aveTCG,
						   double *varPECG, double *varKECG,double *varTCG,
						   struct my_netcdf_out_id_MCD nc_id_MCDCG,  FILE *outputfileCG,
						   // Z  /////////////////////////////////////////////////////////
						   double *Z,double *velZ,double massZ,
						   double *zetaZ,double *V_zetaZ,
						   double QZ,double NfKTZ,double TZ,
						   int numZ_dih,int numZ_ang,int numZ_bon,
						   double KZAA, double KZCG,
						   int **pairs_dihed_AA,int **pairs_dihed_CG,
						   int **pairs_angle_AA,int **pairs_angle_CG,
						   int **pairs_bond_AA, int **pairs_bond_CG,
						   double *avePEZ, double *aveKEZ,double *aveTZ,
						   double *varPEZ, double *varKEZ,double *varTZ, 
						   FILE *trjfileZ, FILE *trjfilThetaAA, FILE *trjfilThetaCG,
						   // CM  /////////////////////////////////////////////////////////
						   double *mass, double *massCA, int numatom, int numCAatom,
						   int numstep, int interval,int *l,
						   double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
						   double UNITT, double k_B,double pi,
						   double *PEZAA, double *PEZCG, double *PEZ,
						   double *ele_ele,    
						   double *ALJ_s,    double *BLJ_s,
						   double *ele_ele_14, 
						   double *ALJ_14_s, double *BLJ_14_s) {
  int i,j,k;

  double PEAA=0.0,KEAA=0.0,EtAA,PEvAA,KEvAA;
  double PECG=0.0,KECG=0.0,EtCG,PEvCG,KEvCG;
  double KEZ,KEvZ,PEvZ,EtZ;

  double *thetaAA,*thetaCG,*frcZ,*fAA,*fCG;
  double summass,COM[3],crd_nc[MAXATOM][3];

  ffLc_calcffandforce_HS(crdAA,numatom,&e,&f,ap_AA,
  			 ele_ele,ALJ_s,BLJ_s,
  			 ele_ele_14,ALJ_14_s,BLJ_14_s);
  GOLM_Clementi_MB_ff_calcff(crdCG,numCAatom,de,d2,&e_CG);

  fAA=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));
  fCG=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));
  frcZ=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));

  thetaAA=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));
  thetaCG=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));

  TACCM_CTheta_Amber_CAGo_MB_2(crdAA,numatom,thetaAA,
			       numZ_dih,pairs_dihed_AA,
			       numZ_ang,pairs_angle_AA,
			       numZ_bon,pairs_bond_AA,
			       pi);
  TACCM_CTheta_Amber_CAGo_MB_2(crdCG,numatom,thetaCG,
			       numZ_dih,pairs_dihed_CG,
			       numZ_ang,pairs_angle_CG,
			       numZ_bon,pairs_bond_CG,
			       pi);

  TACCM_calc_eff_FF_Z_2_Amber_CAGo_MB(Z,numZ_dih,numZ_ang,numZ_bon,
				      thetaAA,KZAA,fAA,PEZAA,UNITT,pi);
  TACCM_calc_eff_FF_Z_2_Amber_CAGo_MB(Z,numZ_dih,numZ_ang,numZ_bon,
				      thetaCG,KZCG,fCG,PEZCG,UNITT,pi);

  for (i=0;i<numZ_dih+numZ_ang+numZ_bon;++i) frcZ[i]=fAA[i]+fCG[i];

  for (i=0;i<numstep;++i) {
    TACCM_CTheta_Amber_CAGo_MB_2(crdAA,numatom,thetaAA,
				 numZ_dih,pairs_dihed_AA,
				 numZ_ang,pairs_angle_AA,
				 numZ_bon,pairs_bond_AA,
				 pi);
    TACCM_CTheta_Amber_CAGo_MB_2(crdCG,numatom,thetaCG,
				 numZ_dih,pairs_dihed_CG,
				 numZ_ang,pairs_angle_CG,
				 numZ_bon,pairs_bond_CG,
				 pi);

    PEAA=TACCM_MD_Propagetor_NH_MP1998_AA_Amber_CAGo_MB_2(crdAA,velAA,mass,zetaAA,V_zetaAA,
							  QAA,NfKTAA,numatom,&KEAA,&KEvAA,&PEvAA,
							  dt,dt2,nc,wdt4,wdt2,
							  &e,&f,ap_AA,Z,numZ_dih,numZ_ang,numZ_bon,
							  thetaAA,KZAA,
							  pairs_dihed_AA,pairs_angle_AA,pairs_bond_AA,
							  PEZAA,pi,
							  ele_ele,ALJ_s,BLJ_s,
							  ele_ele_14,ALJ_14_s,BLJ_14_s);

    PECG=TACCM_MD_Propagetor_NH_MP1998_CG_Amber_CAGo_MB_2(crdCG,velCG,massCA,
    							  zetaCG,V_zetaCG,QCG,
    							  NfKTCG,numCAatom,&KECG,&KEvCG,&PEvCG,
    							  dt,dt2,nc,wdt4,wdt2,
     							  &e_CG,de,d2,
    							  Z,numZ_dih,numZ_ang,numZ_bon,
							  thetaCG,KZCG,
							  pairs_dihed_CG,pairs_angle_CG,pairs_bond_CG,
							  PEZCG,pi);

    TACCM_CGAA_MD_Propagetor_NH_MP1998_Z_Amber_CAGo_MB_2(Z,velZ,massZ,thetaAA,thetaCG,zetaZ,V_zetaZ,
							 QZ,NfKTZ,numZ_dih,numZ_ang,numZ_bon,
							 &KEZ,&KEvZ,&PEvZ,
							 dt,dt2,nc,wdt4,wdt2,KZAA,KZCG,
							 PEZAA,PEZCG,PEZ,frcZ,UNITT,pi);

    if (i%interval==0) {
      KEAA=KEAA/UNITT;     TAA=KEAA/((3*numatom)*k_B)*2.0;
      PEvAA=PEvAA/UNITT;   KEvAA=KEvAA/UNITT;

      KECG=KECG/UNITT;     TCG=KECG/((3*numatom)*k_B)*2.0;
      PEvCG=PEvCG/UNITT;  KEvCG=KEvCG/UNITT;

      PEAA=e.p_d_t+e.p_a_t+e.p_b_t;
      EtAA=PEAA+KEAA+PEvAA+KEvAA;
      fprintf(outputfileAA,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "
	                    ,i+1,PEAA,KEAA,KEvAA,PEvAA,EtAA,TAA);

      PECG=e_CG.p_MB;
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
      myncL_put_crd_MCD(nc_id_MCDAA,*l,crd_nc);

      summass=0.0; for (j=0;j<numCAatom;++j) summass+=massCA[j];
      for (j=0;j<3;++j) COM[j]=0.0;
      for (j=0;j<numCAatom;++j)  for (k=0;k<3;++k) COM[k]+=massCA[j]*crdCG[j*3+k]/summass;
      for (j=0;j<numCAatom;++j)  for (k=0;k<3;++k) crdCG[j*3+k]-=COM[k];
      for (j=0;j<numCAatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crdCG[j*3+k];
      myncL_put_crd_MCD(nc_id_MCDCG,*l,crd_nc);
      ++(*l);

      ///////////////// TACCM //////////////////////
      KEZ=KEZ/UNITT;      TZ=KEZ/((numZ_dih+numZ_ang+numZ_bon)*k_B)*2.0;

      PEvZ=PEvZ/UNITT;      KEvZ=KEvZ/UNITT;

      EtZ=*PEZ+KEZ+PEvZ+KEvZ;
      fprintf(outputfileAA,"%d %e %e %e %e %e %e %e\n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);

      for (j=0;j<numZ_dih+numZ_ang+numZ_bon;++j) fprintf(trjfileZ,"%e ",Z[j]);
      fprintf(trjfileZ,"\n");
      for (j=0;j<numZ_dih+numZ_ang+numZ_bon;++j) fprintf(trjfilThetaAA,"%e ",thetaAA[j]);
      fprintf(trjfilThetaAA,"\n");      
      for (j=0;j<numZ_dih+numZ_ang+numZ_bon;++j) fprintf(trjfilThetaCG,"%e ",thetaCG[j]);
      fprintf(trjfilThetaCG,"\n");      
      ///////////////// TACCM //////////////////////
    }
  }

  return *PEZ;
}

double TACCM_MD_Propagetor_NH_MP1998_AA_Amber_CAGo_MB_2(double *crd,double *vel,double *mass,
							double *zeta,double *V_zeta,double Q,
							double NfKT,int numatom,double *KE,double *KEv,double *PEv,
							double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
							struct potential *e, struct force *f, struct AmberParmL ap,
							double *Z,  int numZ_dih, int numZ_ang, int numZ_bon,
							double *theta,double Kapa,
							int **pairs_dih, int **pairs_ang, int **pairs_bon,
							double *PEZ,double pi,
							double *ele_ele,    double *ALJ_s,    double *BLJ_s,
							double *ele_ele_14, double *ALJ_14_s, double *BLJ_14_s) {
  int i,j,k;
  double *frc;
  double **frcZ;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frcZ=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) frcZ[i]=(double *)gcemalloc(sizeof(double)*3);

  TACCM_CTheta_Amber_CAGo_MB_2(crd,numatom,Z,
  			       numZ_dih,pairs_dih,
  			       numZ_ang,pairs_ang,
  			       numZ_bon,pairs_bon,
  			       pi);

  *PEZ=TACCM_calc_eff_FF_MD_Amber_CAGo_MB(crd,numatom,theta,Z,
  					  numZ_dih,numZ_ang,numZ_bon,
  					  Kapa,frcZ,
  					  pairs_dih,pairs_ang,pairs_bon,pi);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_e[i*3+j]+(*f).f_e_14[i*3+j]+frcZ[i][j];

  MDb_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];

  ffLc_calcffandforce_HS(crd,numatom,e,f,ap,
  			 ele_ele,ALJ_s,BLJ_s,
  			 ele_ele_14,ALJ_14_s,BLJ_14_s);

  TACCM_CTheta_Amber_CAGo_MB_2(crd,numatom,Z,
  			       numZ_dih,pairs_dih,
  			       numZ_ang,pairs_ang,
  			       numZ_bon,pairs_bon,pi);
  *PEZ=TACCM_calc_eff_FF_MD_Amber_CAGo_MB(crd,numatom,theta,Z,
  					  numZ_dih,numZ_ang,numZ_bon,
  					  Kapa,frcZ,
  					  pairs_dih,pairs_ang,pairs_bon,pi);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_e[i*3+j]+(*f).f_e_14[i*3+j]+frcZ[i][j];

  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MDb_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  *KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) *KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return 0.0;
}

double TACCM_MD_Propagetor_NH_MP1998_CG_Amber_CAGo_MB_2(double *crd,double *vel,double *mass,
							double *zeta,double *V_zeta,double Q,
							double NfKT,int numCAatom,
							double *KE,double *KEv,double *PEv,
							double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
							struct potential_GOLM_Clementi_MB *e_CG,
							double de, double d2,
							double *Z, int numZ_dih, int numZ_ang, int numZ_bon,
							double *theta,double Kapa,
							int **pairs_dih, int **pairs_ang, int **pairs_bon, 
							double *PEZ,double pi) {
  int i,j,k;
  double *frc;
  double **frcZ;

  frc=(double *)gcemalloc(sizeof(double)*numCAatom*3);
  frcZ=(double **)gcemalloc(sizeof(double *)*numCAatom);
  for (i=0;i<numCAatom;++i) frcZ[i]=(double *)gcemalloc(sizeof(double)*3);

  //TACCM_CTheta(crd,numCAatom,theta,numZ_dih,pairs_dih,pi);
  TACCM_CTheta_Amber_CAGo_MB_2(crd,numCAatom,Z,
  			       numZ_dih,pairs_dih,
  			       numZ_ang,pairs_ang,
  			       numZ_bon,pairs_bon,pi);
  //*PEZ=TACCM_calc_eff_FF_MD(crd,numCAatom,theta,Z,numZ_dih,Kapa,frcZ,pairs_dih,pi);
  *PEZ=TACCM_calc_eff_FF_MD_Amber_CAGo_MB(crd,numCAatom,theta,Z,
  					  numZ_dih,numZ_ang,numZ_bon,
  					  Kapa,frcZ,
  					  pairs_dih,pairs_ang,pairs_bon,pi);

  for (i=0;i<numCAatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_CG).f_MB[i][j]+frcZ[i][j];

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numCAatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numCAatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numCAatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];

  GOLM_Clementi_MB_ff_calcff(crd,numCAatom,de,d2,e_CG);

  TACCM_CTheta_Amber_CAGo_MB_2(crd,numCAatom,Z,
  			       numZ_dih,pairs_dih,
  			       numZ_ang,pairs_ang,
  			       numZ_bon,pairs_bon,pi);
  //*PEZ=TACCM_calc_eff_FF_MD(crd,numCAatom,theta,Z,numZ_dih,Kapa,frcZ,pairs_dih,pi);
  *PEZ=TACCM_calc_eff_FF_MD_Amber_CAGo_MB(crd,numCAatom,theta,Z,
  					  numZ_dih,numZ_ang,numZ_bon,
  					  Kapa,frcZ,
  					  pairs_dih,pairs_ang,pairs_bon,pi);

  for (i=0;i<numCAatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_CG).f_MB[i][j]+frcZ[i][j];
  for (i=0;i<numCAatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numCAatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numCAatom,nc,wdt4,wdt2);

  *KE=0.0; for (i=0;i<numCAatom;++i) for (j=0;j<3;++j) *KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return 0.0;
}

double TACCM_calc_eff_FF_MD_Amber_CAGo_MB(double *crd,int numatom, double *theta, double *Z,  
					  int numZ_dih,int numZ_ang,int numZ_bon,
					  double Kapa, double **frcZ,
					  int **pairs_dih,int **pairs_ang,int **pairs_bon,
					  double pi){
  int i,j,k,l;
  int ii,jj,kk,ll;

  double f;
  double lenij,lenkj;
  double vij[3],vkj[3];
  double cosijk,angijk;
  double f1,f2;

  double atom[4][3];

  double *dvdpsi;
  double delta;

  double PE=0.0;

  int **pairs_temp;

  /////////////////////////////////////////////////////////////////////////////////////
  dvdpsi=(double *)gcemalloc(sizeof(double)*numZ_dih);
  pairs_temp=(int **)gcemalloc(sizeof(int *)*numZ_dih);
  for (i=0;i<numZ_dih;++i) pairs_temp[i]=(int *)gcemalloc(sizeof(int)*4);

  for (i=0;i<numZ_dih;++i) for (j=0;j<4;++j) pairs_temp[i][j]=pairs_dih[i][j]-1;

  for (i=0;i<numZ_dih;++i) {
    if ((delta=Z[i]-theta[i])>pi) delta-=2.0*pi;
    else if ((delta=Z[i]-theta[i])<-1.0*pi) delta+=2.0*pi;
    dvdpsi[i]=-Kapa*delta;
    PE+=0.5*Kapa*delta*delta;
  }

  FVDIHED_force_dihed(crd,numatom,frcZ,pairs_temp,dvdpsi,numZ_dih);
  /////////////////////////////////////////////////////////////////////////////////////////

  for (i=0;i<numZ_ang;++i) {
    ii=pairs_ang[i][0]-1;
    jj=pairs_ang[i][1]-1;
    kk=pairs_ang[i][2]-1;

    for (j=0;j<3;++j) {
      atom[0][j]=crd[ii*3+j];
      atom[1][j]=crd[jj*3+j];
      atom[2][j]=crd[kk*3+j];
    }

    lenij = len(atom[0],atom[1]);
    lenkj = len(atom[2],atom[1]);
    for (j=0;j<3;++j) {
      vij[j]=atom[1][j]-atom[0][j];
      vkj[j]=atom[1][j]-atom[2][j];
    }
    cosijk=inprod(vij,vkj,3);
    cosijk=cosijk/lenij/lenkj;
    if (cosijk>=1.0) angijk=0.0;
    else if (cosijk<=0.0) angijk=pi;
    else  angijk = acos(cosijk);

    angijk = ang(atom[0],atom[1],atom[2]);

    //    PE += Kapa*(angijk-theta[i+numZ_bon])*(angijk-theta[i+numZ_bon]);

    for (j=0;j<3;++j) {
      f1 = -2.0*Kapa*(angijk-theta[i+numZ_dih])/(lenij*sin(angijk))*(vkj[j]/lenkj-cosijk*vij[j]/lenij)
	*4.184070*100.0;
      f2 = -2.0*Kapa*(angijk-theta[i+numZ_dih])/(lenkj*sin(angijk))*(vij[j]/lenij-cosijk*vkj[j]/lenkj)
	*4.184070*100.0;

      frcZ[ii][j] += f1;
      frcZ[jj][j] += f2;
      frcZ[kk][j] += -f1-f2;
    }
  }
  /////////////////////////////////////////////////////////////////////////////////////////

  for (i=0;i<numZ_bon;++i) {
    ii=pairs_bon[i][0]-1;
    jj=pairs_bon[i][1]-1;

    for (j=0;j<3;++j) {
      atom[0][j]=crd[ii*3+j];
      atom[1][j]=crd[jj*3+j];
    }
  
    lenij = len(atom[0],atom[1]);
    // PE+=Kapa*(lenij-theta[i+numZ_dih+numZ_ang])*(lenij-theta[i+numZ_dih+numZ_ang]);

    for (j=0;j<3;++j) {
      f = 2.0*Kapa*(lenij-theta[i+numZ_dih+numZ_ang])*(atom[1][j]-atom[0][j])/lenij*4.184070*100.0;
      frcZ[ii][j] += f;
      frcZ[jj][j] += -f;
    }
  }

  return PE;
}

double TACCM_CTheta_Amber_CAGo_MB_2(double *crd,int numatom,double *theta, 
				    int numdihe, int **pairs_dih_AA,
				    int numangl, int **pairs_ang_AA,
				    int numbond, int **pairs_bon_AA, 
				    double pi){
  int i,j,k,l;
  int ii,jj,kk,ll;

  double lenij,lenkj;
  double cosijk,angijk;

  double m[3],n[3],m_n[3],n_n[3],lm,ln;
  double vij[3],vkj[3],vkl[3];
  double lkj;
  double vijvkj,vklvkj;

  double atom[4][3];
  double angl,dihed;

  for (i=0;i<numdihe;++i) {
    ii=pairs_dih_AA[i][0]-1;
    jj=pairs_dih_AA[i][1]-1;
    kk=pairs_dih_AA[i][2]-1;
    ll=pairs_dih_AA[i][3]-1;

    for (j=0;j<3;++j) {
      atom[0][j]=crd[ii*3+j];
      atom[1][j]=crd[jj*3+j];
      atom[2][j]=crd[kk*3+j];
      atom[3][j]=crd[ll*3+j];
    }

    for (j=0;j<3;++j) {
      vij[j] = atom[1][j]-atom[0][j];
      vkj[j] = atom[1][j]-atom[2][j];
      vkl[j] = atom[3][j]-atom[2][j];
    }
    lkj=sqrt(inprod(vkj,vkj,3));
    
    outprod(vij,vkj,m);
    outprod(vkj,vkl,n);
    lm=sqrt(inprod(m,m,3));
    ln=sqrt(inprod(n,n,3));
    for (j=0;j<3;++j) {
      m_n[j]=m[j]/lm;
      n_n[j]=n[j]/ln;
    }
    
    dihed=inprod(m_n,n_n,3);
    if (dihed>=1.0)
      dihed=0.0;
    else if (dihed<=-1.0)
      dihed=pi;
    else
      dihed=acos(dihed);
    if (inprod(vij,n,3)>0) dihed=-dihed;
    if (dihed<-1.0*pi) dihed=2.0*pi+dihed;
    if (dihed>pi) dihed=-2.0*pi+dihed;

    theta[i]=dihed;
  }

  ////////////////////////////////////////////////////////////

  for (i=0;i<numangl;++i) {
    ii=pairs_ang_AA[i][0]-1;
    jj=pairs_ang_AA[i][1]-1;
    kk=pairs_ang_AA[i][2]-1;
    for (j=0;j<3;++j) {
      atom[0][j]=crd[ii*3+j];
      atom[1][j]=crd[jj*3+j];
      atom[2][j]=crd[kk*3+j];
    }
  
    lenij = len(atom[0],atom[1]);
    lenkj = len(atom[2],atom[1]);
    for (j=0;j<3;++j) {
      vij[j]=atom[1][j]-atom[0][j];
      vkj[j]=atom[1][j]-atom[2][j];
    }
    cosijk=inprod(vij,vkj,3);
    cosijk=cosijk/lenij/lenkj;
    if (cosijk>=1.0) angijk=0.0;
    else if (cosijk<=0.0) angijk=pi;
    else  angijk = acos(cosijk);
  
    theta[i+numdihe]=angijk;
  }

  ////////////////////////////////////////////////////////////
  for (i=0;i<numbond;++i) {
    ii=pairs_bon_AA[i][0]-1;
    jj=pairs_bon_AA[i][1]-1;
    for (j=0;j<3;++j) {
      atom[0][j]=crd[ii*3+j];
      atom[1][j]=crd[jj*3+j];
    }
    lenij=0.0;
    for (j=0;j<3;++j) {
      lenij += (atom[0][j]-atom[1][j])*(atom[0][j]-atom[1][j]);
    }
    lenij=sqrt(lenij);
    theta[i+numdihe+numangl]=lenij;
  }

  ////////////////////////////////////////////////////////////

  return 0.0;
}

double TACCM_calc_eff_FF_Z_2_Amber_CAGo_MB(double *Z,int numZ_dih, int numZ_ang, int numZ_bon,
					   double *theta,double KZ,double *f, double *PE,
					   double UNITT, double pi){
  int i,j;
  //  double PE=0.0;
  double delta;

  *PE=0.0;
  for (i=0;i<numZ_dih;++i) {
    //    delta=Z[i]-theta[i];
    //    if (delta < 0) delta=-1.0*delta;
    //    if (delta>pi) delta=2.0*pi-delta;
    if ((delta=Z[i]-theta[i])>pi) delta-=2.0*pi;
    else if ((delta=Z[i]-theta[i])<-1.0*pi) delta+=2.0*pi;
    *PE+=0.5*KZ*delta*delta;
    f[i]=-KZ*delta*UNITT;
  }

  for (i=0;i<numZ_ang;++i) {
    //    delta=Z[i]-theta[i];
    //    if (delta < 0) delta=-1.0*delta;
    //    if (delta>pi) delta=2.0*pi-delta;
    if ((delta=Z[i+numZ_dih]-theta[i+numZ_dih])>pi) delta-=2.0*pi;
    else if ((delta=Z[i+numZ_dih]-theta[i+numZ_dih])<-1.0*pi) delta+=2.0*pi;
    *PE+=0.5*KZ*delta*delta;
    f[i+numZ_dih]=-KZ*delta*UNITT;
  }

  for (i=0;i<numZ_bon;++i) {
    delta=Z[i+numZ_dih+numZ_ang]-theta[i+numZ_dih+numZ_ang];
    *PE+=0.5*KZ*delta*delta;
    f[i+numZ_dih+numZ_ang]=-KZ*delta*UNITT;
  }

  return *PE;
}
