
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "FFL.h"

#include "GOLMAA_PROTEINS2008_set.h"
#include "GOLMAA_PROTEINS2008.h"

#include "TACCM_CGAAMDrun_Amber_PROTEINS2008.h"
#include "TACCM_CGFGABAbMDrun_Amber_PROTEINS2008.h"

#include "TACCM_MDrun.h"
#include "TACCM_MD.h"

#include "MD_NHC_MP1996.h"
#include "MDrun.h"
#include "MD.h"

#include "ABAb.h"

double runTACCM_CGFG_ABAbMD_NH_new_Amber_PROTEINS2008(// FG ////////////////////////////////////////////
						      double *FGcrd,double *FGq,double *FGqvel,
						      double *FGpredict,double *FGcorrect,
						      double FGs, double FGs_vel, 
						      double FGpredict_s[5],   double FGcorrect_s[5],
						      double FGgzi, double FGgzi_vel,
						      double FGpredict_gzi[5], double FGcorrect_gzi[5],
						      int *FGnumclutparent,int *FGterminal,int *FGorigin,
						      double *FGvel_Term,
						      double **FGpredict_Term,double **FGpredict_Term2,
						      double **FGcorrect_Term,double **FGcorrect_Term2,
						      CLTb *FGclt, double FGQ,
						      struct potential e, struct force f, 
						      double TFG, double TFGo, double KEFGo,
						      double *avePEFG, double *aveKEFG,double *aveTFG,
						      double *varPEFG, double *varKEFG,double *varTFG,
						      struct my_netcdf_out_id_MCD nc_id_MCDFG,  FILE *outputfileFG,
						      // CG ////////////////////////////////////////////
						      double *CGcrd,double *CGq,double *CGqvel,
						      double *CGpredict,double *CGcorrect,
						      double CGs, double CGs_vel, 
						      double CGpredict_s[5],   double CGcorrect_s[5],
						      double CGgzi, double CGgzi_vel,
						      double CGpredict_gzi[5], double CGcorrect_gzi[5],
						      int *CGnumclutparent,int *CGterminal,int *CGorigin,
						      double *CGvel_Term,
						      double **CGpredict_Term,double **CGpredict_Term2,
						      double **CGcorrect_Term,double **CGcorrect_Term2,
						      CLTb *CGclt, double CGQ,
						      struct potential_GOLMAA_PROTEINS2008 e_CG, 
						      double TCG, double TCGo, double KECGo,
						      double *avePECG, double *aveKECG,double *aveTCG,
						      double *varPECG, double *varKECG,double *varTCG,
						      struct my_netcdf_out_id_MCD nc_id_MCDCG,  FILE *outputfileCG,
						      // Z  /////////////////////////////////////////////
						      double *Z,double *velZ,double massZ,
						      double **predict_Z,double **correct_Z,
						      double sZ,double s_velZ,double gziZ,double gzi_velZ,
						      double predict_gziZ[5],double correct_gziZ[5],
						      double predict_sZ[5],double correct_sZ[5],
						      double QZ,double TZ, double TZo, double KEZo,
						      int numZ,double KZFG, double KZCG,int **pairs,
						      double *avePEZ, double *aveKEZ,double *aveTZ,
						      double *varPEZ, double *varKEZ,double *varTZ, 
						      FILE *trjfileZ, FILE *trjfilThetaFG, FILE *trjfilThetaCG,
						      // CM  ///////////////////////////////////////////////
						      double *mass, int numatom, int numclut,int DOF,
						      int numstep, int interval,int *l,
						      double dt,double tau,double tau2,
						      double UNITT, double k_B,double pi,
						      double *PEZFG, double *PEZCG, double *PEZ) {
  int i,j,k;

  double *Q_frcFCG,*Q_frcFCG_d_Amber,*Q_d_Amber;

  double *Q_frcFG,*frcFG,potFG;
  double /**qFG,*/*qaccFG,*qrotFG;
  double *delta_TermFG,*acc_TermFG,*acc_TermFG2;

  double *Q_frcCG,*frcCG,potCG;
  double /**qCG,*/*qaccCG,*qrotCG;
  double *delta_TermCG,*acc_TermCG,*acc_TermCG2;

  double *accZ;

  double PEFG=0.0,KEFG=0.0,EtFG,PEvFG,KEvFG;
  double PECG=0.0,KECG=0.0,EtCG,PEvCG,KEvCG;
  double KEZ,KEvZ,PEvZ,EtZ;

  double *thetaFG,*thetaCG,*frcZ,*fFG,*fCG,**fFG_MD,**fCG_MD;
  double summass,COM[3],crd_nc[MAXATOM][3];

  FGs=1.0;FGs_vel=0.0;FGgzi=0.0;
  CGs=1.0;CGs_vel=0.0;CGgzi=0.0;
  sZ=1.0;s_velZ=0.0;gziZ=0.0;

  /****************************************************/
  /* qFG=(double *)gcemalloc(sizeof(double)*numclut); */
  /* qCG=(double *)gcemalloc(sizeof(double)*numclut); */
  /* for (i=0;i<numclut;++i) {			      */
  /*   qFG[i]=0.0;				      */
  /*   qCG[i]=0.0;				      */
  /* }						      */
  /****************************************************/

  ABAbNH_set_new(FGs,FGs_vel,FGgzi,FGpredict_gzi,FGcorrect_gzi,FGpredict_s,FGcorrect_s,tau,&tau2,&FGQ,KEFGo,dt);
  ABAb_integ_set(FGq,FGqvel,FGpredict,FGcorrect,numclut,dt);

  ABAbNH_set_new(CGs,CGs_vel,CGgzi,CGpredict_gzi,CGcorrect_gzi,CGpredict_s,CGcorrect_s,tau,&tau2,&CGQ,KECGo,dt);
  ABAb_integ_set(CGq,CGqvel,CGpredict,CGcorrect,numclut,dt);

  TACCM_NH_set_new(sZ,s_velZ,gziZ,predict_gziZ,correct_gziZ,predict_sZ,correct_sZ,tau,&tau2,&QZ,KEZo,dt);

  Q_frcFCG=(double *)gcemalloc(sizeof(double)*numclut);
  Q_frcFCG_d_Amber=(double *)gcemalloc(sizeof(double)*numclut);
  Q_d_Amber=(double *)gcemalloc(sizeof(double)*numclut);

  Q_frcFG=(double *)gcemalloc(sizeof(double)*numclut);
  frcFG=(double *)gcemalloc(sizeof(double)*numatom*3);
  qrotFG=(double *)gcemalloc(sizeof(double)*numclut);
  qaccFG=(double *)gcemalloc(sizeof(double)*numclut);
  delta_TermFG=(double *)gcemalloc(sizeof(double)*6);
  acc_TermFG=(double *)gcemalloc(sizeof(double)*6);
  acc_TermFG2=(double *)gcemalloc(sizeof(double)*6);

  Q_frcCG=(double *)gcemalloc(sizeof(double)*numclut);
  frcCG=(double *)gcemalloc(sizeof(double)*numatom*3);
  qrotCG=(double *)gcemalloc(sizeof(double)*numclut);
  qaccCG=(double *)gcemalloc(sizeof(double)*numclut);
  delta_TermCG=(double *)gcemalloc(sizeof(double)*6);
  acc_TermCG=(double *)gcemalloc(sizeof(double)*6);
  acc_TermCG2=(double *)gcemalloc(sizeof(double)*6);

  ffL_calcffandforce(FGcrd,numatom,&e,&f);
  GOLMAA_PROTEINS2008_ff_calcff_wobaimp_b(CGcrd,numatom,&e_CG);

  fFG=(double *)gcemalloc(sizeof(double)*numZ);
  fCG=(double *)gcemalloc(sizeof(double)*numZ);
  frcZ=(double *)gcemalloc(sizeof(double)*numZ);
  accZ=(double *)gcemalloc(sizeof(double)*numZ);

  fFG_MD=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) fFG_MD[i]=(double *)gcemalloc(sizeof(double)*3);
  fCG_MD=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) fCG_MD[i]=(double *)gcemalloc(sizeof(double)*3);

  thetaFG=(double *)gcemalloc(sizeof(double)*numZ);
  thetaCG=(double *)gcemalloc(sizeof(double)*numZ);

  TACCM_CTheta(FGcrd,numatom,thetaFG,numZ,pairs,pi);
  TACCM_CTheta(CGcrd,numatom,thetaCG,numZ,pairs,pi);
  TACCM_calc_eff_FF_Z(Z,numZ,thetaFG,KZFG,fFG,pi);
  TACCM_calc_eff_FF_Z(Z,numZ,thetaCG,KZCG,fCG,pi);

  for (i=0;i<numZ;++i) frcZ[i]=fFG[i]+fCG[i];

  ABAbNH_calcKE_new(FGgzi,FGs,FGs_vel,FGQ,KEFGo,&PEvFG,&KEvFG);
  ABAbNH_calcKE_new(CGgzi,CGs,CGs_vel,CGQ,KECGo,&PEvCG,&KEvCG);

  for (i=0;i<numstep;++i) {
    ABAb_integ_pret(qrotFG,FGqvel,FGq,FGpredict,FGcorrect,dt,numclut);
    ABAb_integ_pret_Term(FGpredict_Term,FGpredict_Term2,FGcorrect_Term,FGcorrect_Term2,FGvel_Term,delta_TermFG,dt);
    ABAbNH_update_pret_new(&FGgzi,FGpredict_gzi,FGcorrect_gzi,&FGs,&FGs_vel,FGpredict_s,FGcorrect_s,dt);

    ABAb_integ_pret(qrotCG,CGqvel,CGq,CGpredict,CGcorrect,dt,numclut);
    ABAb_integ_pret_Term(CGpredict_Term,CGpredict_Term2,CGcorrect_Term,CGcorrect_Term2,CGvel_Term,delta_TermCG,dt);
    ABAbNH_update_pret_new(&CGgzi,CGpredict_gzi,CGcorrect_gzi,&CGs,&CGs_vel,CGpredict_s,FGcorrect_s,dt);
    
    TACCM_integ_pret_Z(predict_Z,correct_Z,Z,velZ,numZ,dt,pi);
    TACCM_NH_update_pret_new(&gziZ,&gziZ,predict_gziZ,correct_gziZ,&sZ,&s_velZ,predict_sZ,correct_sZ,dt);

    ABAb_update(FGclt,FGcrd,qrotFG,numclut,numatom);
    ABAb_update(CGclt,CGcrd,qrotCG,numclut,numatom);

    for (j=0;j<numclut;++j) Q_frcFG[j]=0.0;
    for (j=0;j<numatom;++j) for (k=0;k<3;++k) frcFG[j*3+k]=0.0;

    for (j=0;j<numclut;++j) Q_frcCG[j]=0.0;
    for (j=0;j<numatom;++j) for (k=0;k<3;++k) frcCG[j*3+k]=0.0;

    for (j=0;j<numclut;++j) Q_d_Amber[j]=0.0;
    for (j=0;j<numclut;++j) Q_frcFCG_d_Amber[j]=0.0;

    ffL_calcTorque(Q_d_Amber,FGcrd,numclut,FGnumclutparent,FGterminal,FGorigin);
    ffL_calcffandforce(FGcrd,numatom,&e,&f);
    for (j=0;j<numatom;++j)for (k=0;k<3;++k)frcFG[j*3+k]=f.f_e[j*3+k]+f.f_LJ[j*3+k]+f.f_e_14[j*3+k]+f.f_LJ_14[j*3+k];

    GOLMAA_PROTEINS2008_ff_calcff_wobaimp_b(CGcrd,numatom,&e_CG);
    for (j=0;j<numatom;++j)for (k=0;k<3;++k) frcCG[j*3+k]=(e_CG).f_t[j][k];

    TACCM_CTheta(FGcrd,numatom,thetaFG,numZ,pairs,pi);
    TACCM_CTheta(CGcrd,numatom,thetaCG,numZ,pairs,pi);
    
    //    *PEZFG=TACCM_calc_eff_FF_MD(FGcrd,numatom,thetaFG,Z,numZ,KZFG,fFG_MD,pairs,pi);
    //    *PEZCG=TACCM_calc_eff_FF_MD(CGcrd,numatom,thetaCG,Z,numZ,KZCG,fCG_MD,pairs,pi);

    //    for (j=0;j<numatom;++j)for (k=0;k<3;++k) frcFG[j*3+k]+=fFG_MD[j][k];
    //    for (j=0;j<numatom;++j)for (k=0;k<3;++k) frcCG[j*3+k]+=fCG_MD[j][k];
    
    for (j=0;j<numZ;++j) {
      while (Z[j]>pi) Z[j]-=2.0*pi;
      while (Z[j]<=-pi) Z[j]+=2.0*pi;
    }
    TACCM_calc_eff_FF_Z_2(Z,numZ,thetaFG,KZFG,fFG,PEZFG,pi);
    TACCM_calc_eff_FF_Z_2(Z,numZ,thetaCG,KZCG,fCG,PEZCG,pi);

    /**PEZFG=*/TACCM_calc_eff_FF(thetaFG,Z,numZ,KZFG,Q_frcFG,pairs,pi);
    /**PEZCG=*/TACCM_calc_eff_FF(thetaCG,Z,numZ,KZCG,Q_frcCG,pairs,pi);
    
    //    /**PEZFG=*/TACCM_calc_eff_FF_Z(Z,numZ,thetaFG,KZFG,fFG,pi);
    //    /**PEZCG=*/TACCM_calc_eff_FF_Z(Z,numZ,thetaCG,KZCG,fCG,pi);
    
    for (j=0;j<numclut;++j) Q_frcFCG[j]=0.0;

    *PEZ=(*PEZFG)+(*PEZCG);
    for (j=0;j<numZ;++j) frcZ[j]=fFG[j]+fCG[j];
    for (j=0;j<numclut;++j) Q_frcFCG_d_Amber[j]=Q_frcFG[j]+Q_d_Amber[j];

    ABAb_calcKineE_TermOn(&KEFG,&KEvFG,&PEvFG,KEFGo,FGclt,FGcrd,FGqvel,FGs,FGs_vel,
			  FGQ,FGvel_Term,numclut,numatom,NVT,numclut+6);

    ABAb_calcKineE_TermOn(&KECG,&KEvCG,&PEvCG,KECGo,CGclt,CGcrd,CGqvel,CGs,CGs_vel,
			  CGQ,CGvel_Term,numclut,numatom,NVT,numclut+6);

    TFG=KEFG/(DOF*k_B)*2.0;
    TCG=KECG/(DOF*k_B)*2.0;

    TACCM_calcKineE_Z(&KEZ,massZ,velZ,numZ);
    TZ=KEZ/(numZ*k_B)*2.0;

    solverABAb_TermOn_NH_new(qaccFG,FGqvel,FGclt,Q_frcFCG_d_Amber,frcFG,FGcrd,
			     numclut,numatom,FGq,FGgzi,&(FGgzi_vel),FGs,&(FGs_vel),
			     tau2,acc_TermFG,acc_TermFG2,FGvel_Term,TFG,TFGo);

    solverABAb_TermOn_NH_new(qaccCG,CGqvel,CGclt,Q_frcCG,frcCG,CGcrd,
			     numclut,numatom,CGq,CGgzi,&(CGgzi_vel),CGs,&(CGs_vel),
			     tau2,acc_TermCG,acc_TermCG2,CGvel_Term,TCG,TCGo);

    TACCM_solver_NH_Z(accZ,velZ,massZ,frcZ,numZ,gziZ,&gzi_velZ,tau2,TZ,TZo);

    ABAb_integ_cort(qrotFG,FGqvel,FGq,qaccFG,FGpredict,FGcorrect,dt,numclut);
    ABAb_integ_cort_Term(FGpredict_Term,FGpredict_Term2,
			 FGcorrect_Term,FGcorrect_Term2,
			 acc_TermFG,acc_TermFG2,FGvel_Term,
			 delta_TermFG,dt);
    ABAbNH_update_cort_new(&FGgzi,FGgzi_vel,&FGs,&FGs_vel,
			   FGpredict_gzi,FGcorrect_gzi,
			   FGpredict_s,FGcorrect_s,dt);
    ABAb_update(FGclt,FGcrd,qrotFG,numclut,numatom);

    ABAb_integ_cort(qrotCG,CGqvel,CGq,qaccCG,CGpredict,CGcorrect,dt,numclut);
    ABAb_integ_cort_Term(CGpredict_Term,CGpredict_Term2,
			 CGcorrect_Term,CGcorrect_Term2,
			 acc_TermCG,acc_TermCG2,CGvel_Term,
			 delta_TermCG,dt);
    ABAbNH_update_cort_new(&CGgzi,CGgzi_vel,&CGs,&CGs_vel,
			   CGpredict_gzi,CGcorrect_gzi,
			   CGpredict_s,CGcorrect_s,dt);
    ABAb_update(CGclt,CGcrd,qrotCG,numclut,numatom);

    TACCM_integ_cort_Z(predict_Z,correct_Z,accZ,Z,velZ,numZ,dt,pi);
    TACCM_NH_update_cort_new(&gziZ,gzi_velZ,&sZ,&s_velZ,predict_gziZ,correct_gziZ,predict_sZ,correct_sZ,dt);

    ABAb_calcKineE_TermOn_new_simp(&KEFG,FGclt,FGcrd,FGqvel,FGvel_Term,numclut,numatom);
    ABAbNH_calcKE_new(FGgzi,FGs,FGs_vel,FGQ,KEFGo,&PEvFG,&KEvFG);

    ABAb_calcKineE_TermOn_new_simp(&KECG,CGclt,CGcrd,CGqvel,CGvel_Term,numclut,numatom);
    ABAbNH_calcKE_new(CGgzi,CGs,CGs_vel,CGQ,KECGo,&PEvCG,&KEvCG);

    TACCM_calcKineE_Z(&KEZ,massZ,velZ,numZ);
    TACCM_NH_calcKE_new(gziZ,sZ,s_velZ,/*Q_NH*/QZ,/*KEobj*/KEZo,&PEvZ,&KEvZ);
    TZ=KEZ/(numZ*k_B)*2.0;

    if (i%interval==0) {
      //      KEFG=KEFG/UNITT;     TFG=KEFG/(DOF*k_B)*2.0;
      //      PEvFG=PEvFG/UNITT;   KEvFG=KEvFG/UNITT;
      
      //      KECG=KECG/UNITT;     TCG=KECG/(DOF*k_B)*2.0;
      //      PEvCG=PEvCG/UNITT;  KEvCG=KEvCG/UNITT;

      PEFG=0.5*e.p_e_t+0.5*e.p_LJ_t+0.5*e.p_e_14_t+0.5*e.p_LJ_14_t+e.p_d_t/*+e.p_a_t+e.p_b_t*/;
      EtFG=PEFG+KEFG+PEvFG+KEvFG;
      fprintf(outputfileFG,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "
                  	   ,i+1,PEFG, KEFG, KEvFG,PEvFG,EtFG,*PEZFG,TFG);

      PECG=e_CG.p_t;
      EtCG=PECG+KECG+PEvCG+KEvCG;
      fprintf(outputfileCG,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "
  	                   ,i+1,PECG,KECG, KEvCG, PEvCG,EtCG,*PEZCG,TCG);

      summass=0.0; for (j=0;j<numatom;++j) summass+=mass[j];
      for (j=0;j<3;++j) COM[j]=0.0;
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) COM[k]+=mass[j]*FGcrd[j*3+k]/summass;
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) FGcrd[j*3+k]-=COM[k];
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=FGcrd[j*3+k];
      myncL_put_crd_ene_MCD(nc_id_MCDFG,*l,crd_nc,e,0.0);

      summass=0.0; for (j=0;j<numatom;++j) summass+=mass[j];
      for (j=0;j<3;++j) COM[j]=0.0;
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) COM[k]+=mass[j]*CGcrd[j*3+k]/summass;
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) CGcrd[j*3+k]-=COM[k];
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=CGcrd[j*3+k];
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=CGcrd[j*3+k];
      myncL_put_crd_ene_MCD(nc_id_MCDCG,*l,crd_nc,e,0.0);
      ++(*l);

      ///////////////// TACCM //////////////////////
      //      KEZ=KEZ/UNITT;      TZ=KEZ/(numZ*k_B)*2.0;

      //      PEvZ=PEvZ/UNITT;      KEvZ=KEvZ/UNITT;

      EtZ=*PEZ+KEZ+PEvZ+KEvZ;
      fprintf(outputfileFG,"%d %e %e %e %e %e %e \n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);

      fprintf(outputfileCG,"%d %e %e %e %e %e %e \n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);

      for (j=0;j<numZ;++j) fprintf(trjfileZ,"%e ",Z[j]);
      fprintf(trjfileZ,"\n");
      for (j=0;j<numZ;++j) fprintf(trjfilThetaFG,"%e ",thetaFG[j]);
      fprintf(trjfilThetaFG,"\n");      
      for (j=0;j<numZ;++j) fprintf(trjfilThetaCG,"%e ",thetaCG[j]);
      fprintf(trjfilThetaCG,"\n");      
      ///////////////// TACCM //////////////////////
    }
  }

  return *PEZ;
}
