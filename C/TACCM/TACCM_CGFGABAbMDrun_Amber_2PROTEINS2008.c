
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
#include "TACCM_CGFGABAbMDrun_Amber_2PROTEINS2008.h"

#include "TACCM_MDrun.h"
#include "TACCM_MD.h"

#include "MD_NHC_MP1996.h"
#include "MDrun.h"
#include "MD.h"

#include "ABAb.h"

double runTACCM_2CG1FG_ABAbMD_NH_new_Amber_PROTEINS2008(// FG ////////////////////////////////////////////
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
							// CG1 ////////////////////////////////////////////
							double *CG1crd,double *CG1q,double *CG1qvel,
							double *CG1predict,double *CG1correct,
							double CG1s, double CG1s_vel, 
							double CG1predict_s[5],   double CG1correct_s[5],
							double CG1gzi, double CG1gzi_vel,
							double CG1predict_gzi[5], double CG1correct_gzi[5],
							int *CG1numclutparent,int *CG1terminal,int *CG1origin,
							double *CG1vel_Term,
							double **CG1predict_Term,double **CG1predict_Term2,
							double **CG1correct_Term,double **CG1correct_Term2,
							CLTb *CG1clt, double CG1Q,
							struct potential_GOLMAA_PROTEINS2008 e_CG1, 
							double TCG1, double TCG1o, double KECG1o,
							double *avePECG1, double *aveKECG1,double *aveTCG1,
							double *varPECG1, double *varKECG1,double *varTCG1,
							struct my_netcdf_out_id_MCD nc_id_MCDCG1,  FILE *outputfileCG1,
							// CG2 ////////////////////////////////////////////
							double *CG2crd,double *CG2q,double *CG2qvel,
							double *CG2predict,double *CG2correct,
							double CG2s, double CG2s_vel, 
							double CG2predict_s[5],   double CG2correct_s[5],
							double CG2gzi, double CG2gzi_vel,
							double CG2predict_gzi[5], double CG2correct_gzi[5],
							int *CG2numclutparent,int *CG2terminal,int *CG2origin,
							double *CG2vel_Term,
							double **CG2predict_Term,double **CG2predict_Term2,
							double **CG2correct_Term,double **CG2correct_Term2,
							CLTb *CG2clt, double CG2Q,
							struct potential_GOLMAA_PROTEINS2008 e_CG2, 
							double TCG2, double TCG2o, double KECG2o,
							double *avePECG2, double *aveKECG2,double *aveTCG2,
							double *varPECG2, double *varKECG2,double *varTCG2,
							struct my_netcdf_out_id_MCD nc_id_MCDCG2,  FILE *outputfileCG2,
							// Z  /////////////////////////////////////////////
							double *Z,double *velZ,double massZ,
							double **predict_Z,double **correct_Z,
							double sZ,double s_velZ,double gziZ,double gzi_velZ,
							double predict_gziZ[5],double correct_gziZ[5],
							double predict_sZ[5],double correct_sZ[5],
							double QZ,double TZ, double TZo, double KEZo,
							int numZ,double KZFG, 
							double KZCG1, double KZCG2,int **pairs,
							double *avePEZ, double *aveKEZ,double *aveTZ,
							double *varPEZ, double *varKEZ,double *varTZ, 
							FILE *trjfileZ, FILE *trjfilThetaFG, 
							FILE *trjfilThetaCG1,FILE *trjfilThetaCG2,
							// CM  ///////////////////////////////////////////////
							double *mass, int numatom, int numclut,int DOF,
							int numstep, int interval,int *l,
							double dt,double tau,double tau2,
							double UNITT, double k_B,double pi,
							double *PEZFG, 
							double *PEZCG1, double *PEZCG2, double *PEZ) {
  int i,j,k;

  double *Q_frcFCG,*Q_frcFCG_d_Amber,*Q_d_Amber;

  double *Q_frcFG,*frcFG,potFG;
  double /**qFG,*/*qaccFG,*qrotFG;
  double *delta_TermFG,*acc_TermFG,*acc_TermFG2;

  double *Q_frcCG1,*frcCG1,potCG1;
  double /**qCG1,*/*qaccCG1,*qrotCG1;
  double *delta_TermCG1,*acc_TermCG1,*acc_TermCG12;

  double *Q_frcCG2,*frcCG2,potCG2;
  double /**qCG2,*/*qaccCG2,*qrotCG2;
  double *delta_TermCG2,*acc_TermCG2,*acc_TermCG22;

  double *accZ;

  double PEFG=0.0,KEFG=0.0,EtFG,PEvFG,KEvFG;
  double PECG1=0.0,KECG1=0.0,EtCG1,PEvCG1,KEvCG1;
  double PECG2=0.0,KECG2=0.0,EtCG2,PEvCG2,KEvCG2;
  double KEZ,KEvZ,PEvZ,EtZ;

  double *thetaFG,*thetaCG1,*thetaCG2,*frcZ,*fFG,*fCG1,*fCG2,**fFG_MD,**fCG1_MD,**fCG2_MD;
  double summass,COM[3],crd_nc[MAXATOM][3];

  FGs=1.0;FGs_vel=0.0;FGgzi=0.0;
  CG1s=1.0;CG1s_vel=0.0;CG1gzi=0.0;
  CG2s=1.0;CG2s_vel=0.0;CG2gzi=0.0;
  sZ=1.0;s_velZ=0.0;gziZ=0.0;

  /****************************************************/
  /* qFG=(double *)gcemalloc(sizeof(double)*numclut); */
  /* qCG1=(double *)gcemalloc(sizeof(double)*numclut); */
  /* for (i=0;i<numclut;++i) {			      */
  /*   qFG[i]=0.0;				      */
  /*   qCG1[i]=0.0;				      */
  /* }						      */
  /****************************************************/

  ABAbNH_set_new(FGs,FGs_vel,FGgzi,FGpredict_gzi,FGcorrect_gzi,FGpredict_s,FGcorrect_s,tau,&tau2,&FGQ,KEFGo,dt);
  ABAb_integ_set(FGq,FGqvel,FGpredict,FGcorrect,numclut,dt);

  ABAbNH_set_new(CG1s,CG1s_vel,CG1gzi,CG1predict_gzi,CG1correct_gzi,CG1predict_s,CG1correct_s,tau,&tau2,&CG1Q,KECG1o,dt);
  ABAb_integ_set(CG1q,CG1qvel,CG1predict,CG1correct,numclut,dt);

  ABAbNH_set_new(CG2s,CG2s_vel,CG2gzi,CG2predict_gzi,CG2correct_gzi,CG2predict_s,CG2correct_s,tau,&tau2,&CG2Q,KECG2o,dt);
  ABAb_integ_set(CG2q,CG2qvel,CG2predict,CG2correct,numclut,dt);

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

  Q_frcCG1=(double *)gcemalloc(sizeof(double)*numclut);
  frcCG1=(double *)gcemalloc(sizeof(double)*numatom*3);
  qrotCG1=(double *)gcemalloc(sizeof(double)*numclut);
  qaccCG1=(double *)gcemalloc(sizeof(double)*numclut);
  delta_TermCG1=(double *)gcemalloc(sizeof(double)*6);
  acc_TermCG1=(double *)gcemalloc(sizeof(double)*6);
  acc_TermCG12=(double *)gcemalloc(sizeof(double)*6);

  Q_frcCG2=(double *)gcemalloc(sizeof(double)*numclut);
  frcCG2=(double *)gcemalloc(sizeof(double)*numatom*3);
  qrotCG2=(double *)gcemalloc(sizeof(double)*numclut);
  qaccCG2=(double *)gcemalloc(sizeof(double)*numclut);
  delta_TermCG2=(double *)gcemalloc(sizeof(double)*6);
  acc_TermCG2=(double *)gcemalloc(sizeof(double)*6);
  acc_TermCG22=(double *)gcemalloc(sizeof(double)*6);

  ffL_calcffandforce(FGcrd,numatom,&e,&f);
  GOLMAA_PROTEINS2008_ff_calcff_wobaimp_b(CG1crd,numatom,&e_CG1);
  GOLMAA_PROTEINS2008_ff_calcff_wobaimp_b(CG2crd,numatom,&e_CG2);

  fFG=(double *)gcemalloc(sizeof(double)*numZ);
  fCG1=(double *)gcemalloc(sizeof(double)*numZ);
  fCG2=(double *)gcemalloc(sizeof(double)*numZ);
  frcZ=(double *)gcemalloc(sizeof(double)*numZ);
  accZ=(double *)gcemalloc(sizeof(double)*numZ);

  fFG_MD=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) fFG_MD[i]=(double *)gcemalloc(sizeof(double)*3);
  fCG1_MD=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) fCG1_MD[i]=(double *)gcemalloc(sizeof(double)*3);
  fCG2_MD=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) fCG2_MD[i]=(double *)gcemalloc(sizeof(double)*3);

  thetaFG=(double *)gcemalloc(sizeof(double)*numZ);
  thetaCG1=(double *)gcemalloc(sizeof(double)*numZ);
  thetaCG2=(double *)gcemalloc(sizeof(double)*numZ);

  TACCM_CTheta(FGcrd,numatom,thetaFG,numZ,pairs,pi);
  TACCM_CTheta(CG1crd,numatom,thetaCG1,numZ,pairs,pi);
  TACCM_CTheta(CG2crd,numatom,thetaCG2,numZ,pairs,pi);
  TACCM_calc_eff_FF_Z(Z,numZ,thetaFG,KZFG,fFG,pi);
  TACCM_calc_eff_FF_Z(Z,numZ,thetaCG1,KZCG1,fCG1,pi);
  TACCM_calc_eff_FF_Z(Z,numZ,thetaCG2,KZCG2,fCG2,pi);

  for (i=0;i<numZ;++i) frcZ[i]=fFG[i]+fCG1[i]+fCG2[i];

  ABAbNH_calcKE_new(FGgzi,FGs,FGs_vel,FGQ,KEFGo,&PEvFG,&KEvFG);
  ABAbNH_calcKE_new(CG1gzi,CG1s,CG1s_vel,CG1Q,KECG1o,&PEvCG1,&KEvCG1);
  ABAbNH_calcKE_new(CG2gzi,CG2s,CG2s_vel,CG2Q,KECG2o,&PEvCG2,&KEvCG2);

  for (i=0;i<numstep;++i) {
    ABAb_integ_pret(qrotFG,FGqvel,FGq,FGpredict,FGcorrect,dt,numclut);
    ABAb_integ_pret_Term(FGpredict_Term,FGpredict_Term2,FGcorrect_Term,FGcorrect_Term2,FGvel_Term,delta_TermFG,dt);
    ABAbNH_update_pret_new(&FGgzi,FGpredict_gzi,FGcorrect_gzi,&FGs,&FGs_vel,FGpredict_s,FGcorrect_s,dt);

    ABAb_integ_pret(qrotCG1,CG1qvel,CG1q,CG1predict,CG1correct,dt,numclut);
    ABAb_integ_pret_Term(CG1predict_Term,CG1predict_Term2,CG1correct_Term,CG1correct_Term2,CG1vel_Term,delta_TermCG1,dt);
    ABAbNH_update_pret_new(&CG1gzi,CG1predict_gzi,CG1correct_gzi,&CG1s,&CG1s_vel,CG1predict_s,FGcorrect_s,dt);

    ABAb_integ_pret(qrotCG2,CG2qvel,CG2q,CG2predict,CG2correct,dt,numclut);
    ABAb_integ_pret_Term(CG2predict_Term,CG2predict_Term2,CG2correct_Term,CG2correct_Term2,CG2vel_Term,delta_TermCG2,dt);
    ABAbNH_update_pret_new(&CG2gzi,CG2predict_gzi,CG2correct_gzi,&CG2s,&CG2s_vel,CG2predict_s,FGcorrect_s,dt);
    
    TACCM_integ_pret_Z(predict_Z,correct_Z,Z,velZ,numZ,dt,pi);
    TACCM_NH_update_pret_new(&gziZ,&gziZ,predict_gziZ,correct_gziZ,&sZ,&s_velZ,predict_sZ,correct_sZ,dt);

    ABAb_update(FGclt,FGcrd,qrotFG,numclut,numatom);
    ABAb_update(CG1clt,CG1crd,qrotCG1,numclut,numatom);
    ABAb_update(CG2clt,CG2crd,qrotCG2,numclut,numatom);

    for (j=0;j<numclut;++j) Q_frcFG[j]=0.0;
    for (j=0;j<numatom;++j) for (k=0;k<3;++k) frcFG[j*3+k]=0.0;

    for (j=0;j<numclut;++j) Q_frcCG1[j]=0.0;
    for (j=0;j<numatom;++j) for (k=0;k<3;++k) frcCG1[j*3+k]=0.0;

    for (j=0;j<numclut;++j) Q_frcCG2[j]=0.0;
    for (j=0;j<numatom;++j) for (k=0;k<3;++k) frcCG2[j*3+k]=0.0;

    for (j=0;j<numclut;++j) Q_d_Amber[j]=0.0;
    for (j=0;j<numclut;++j) Q_frcFCG_d_Amber[j]=0.0;

    ffL_calcTorque(Q_d_Amber,FGcrd,numclut,FGnumclutparent,FGterminal,FGorigin);
    ffL_calcffandforce(FGcrd,numatom,&e,&f);
    for (j=0;j<numatom;++j)for (k=0;k<3;++k)frcFG[j*3+k]=f.f_e[j*3+k]+f.f_LJ[j*3+k]+f.f_e_14[j*3+k]+f.f_LJ_14[j*3+k];

    GOLMAA_PROTEINS2008_ff_calcff_wobaimp_b(CG1crd,numatom,&e_CG1);
    for (j=0;j<numatom;++j)for (k=0;k<3;++k) frcCG1[j*3+k]=(e_CG1).f_t[j][k];
    //    for (j=0;j<numatom;++j)for (k=0;k<3;++k) frcCG1[j*3+k]=0.0;

    GOLMAA_PROTEINS2008_ff_calcff_wobaimp_b(CG2crd,numatom,&e_CG2);
    for (j=0;j<numatom;++j)for (k=0;k<3;++k) frcCG2[j*3+k]=(e_CG2).f_t[j][k];

    TACCM_CTheta(FGcrd,numatom,thetaFG,numZ,pairs,pi);
    TACCM_CTheta(CG1crd,numatom,thetaCG1,numZ,pairs,pi);
    TACCM_CTheta(CG2crd,numatom,thetaCG2,numZ,pairs,pi);
    
    //    *PEZFG=TACCM_calc_eff_FF_MD(FGcrd,numatom,thetaFG,Z,numZ,KZFG,fFG_MD,pairs,pi);
    //    *PEZCG1=TACCM_calc_eff_FF_MD(CG1crd,numatom,thetaCG1,Z,numZ,KZCG1,fCG1_MD,pairs,pi);

    //    for (j=0;j<numatom;++j)for (k=0;k<3;++k) frcFG[j*3+k]+=fFG_MD[j][k];
    //    for (j=0;j<numatom;++j)for (k=0;k<3;++k) frcCG1[j*3+k]+=fCG1_MD[j][k];
    
    for (j=0;j<numZ;++j) {
      while (Z[j]>pi) Z[j]-=2.0*pi;
      while (Z[j]<=-pi) Z[j]+=2.0*pi;
    }
    TACCM_calc_eff_FF_Z_2(Z,numZ,thetaFG,KZFG,fFG,PEZFG,pi);
    TACCM_calc_eff_FF_Z_2(Z,numZ,thetaCG1,KZCG1,fCG1,PEZCG1,pi);
    TACCM_calc_eff_FF_Z_2(Z,numZ,thetaCG2,KZCG2,fCG2,PEZCG2,pi);

    /**PEZFG=*/TACCM_calc_eff_FF(thetaFG,Z,numZ,KZFG,Q_frcFG,pairs,pi);
    /**PEZCG1=*/TACCM_calc_eff_FF(thetaCG1,Z,numZ,KZCG1,Q_frcCG1,pairs,pi);
    /**PEZCG2=*/TACCM_calc_eff_FF(thetaCG2,Z,numZ,KZCG2,Q_frcCG2,pairs,pi);
    
    //    /**PEZFG=*/TACCM_calc_eff_FF_Z(Z,numZ,thetaFG,KZFG,fFG,pi);
    //    /**PEZCG1=*/TACCM_calc_eff_FF_Z(Z,numZ,thetaCG1,KZCG1,fCG1,pi);
    
    for (j=0;j<numclut;++j) Q_frcFCG[j]=0.0;

    *PEZ=(*PEZFG)+(*PEZCG1)+(*PEZCG2);
    for (j=0;j<numZ;++j) frcZ[j]=fFG[j]+fCG1[j]+fCG2[j];
    for (j=0;j<numclut;++j) Q_frcFCG_d_Amber[j]=Q_frcFG[j]+Q_d_Amber[j];

    ABAb_calcKineE_TermOn(&KEFG,&KEvFG,&PEvFG,KEFGo,FGclt,FGcrd,FGqvel,FGs,FGs_vel,
			  FGQ,FGvel_Term,numclut,numatom,NVT,numclut+6);

    ABAb_calcKineE_TermOn(&KECG1,&KEvCG1,&PEvCG1,KECG1o,CG1clt,CG1crd,CG1qvel,CG1s,CG1s_vel,
			  CG1Q,CG1vel_Term,numclut,numatom,NVT,numclut+6);

    ABAb_calcKineE_TermOn(&KECG2,&KEvCG2,&PEvCG2,KECG2o,CG2clt,CG2crd,CG2qvel,CG2s,CG2s_vel,
			  CG2Q,CG2vel_Term,numclut,numatom,NVT,numclut+6);

    TFG=KEFG/(DOF*k_B)*2.0;
    TCG1=KECG1/(DOF*k_B)*2.0;
    TCG2=KECG2/(DOF*k_B)*2.0;

    TACCM_calcKineE_Z(&KEZ,massZ,velZ,numZ);
    TZ=KEZ/(numZ*k_B)*2.0;

    solverABAb_TermOn_NH_new(qaccFG,FGqvel,FGclt,Q_frcFCG_d_Amber,frcFG,FGcrd,
			     numclut,numatom,FGq,FGgzi,&(FGgzi_vel),FGs,&(FGs_vel),
			     tau2,acc_TermFG,acc_TermFG2,FGvel_Term,TFG,TFGo);

    solverABAb_TermOn_NH_new(qaccCG1,CG1qvel,CG1clt,Q_frcCG1,frcCG1,CG1crd,
			     numclut,numatom,CG1q,CG1gzi,&(CG1gzi_vel),CG1s,&(CG1s_vel),
			     tau2,acc_TermCG1,acc_TermCG12,CG1vel_Term,TCG1,TCG1o);

    solverABAb_TermOn_NH_new(qaccCG2,CG2qvel,CG2clt,Q_frcCG2,frcCG2,CG2crd,
			     numclut,numatom,CG2q,CG2gzi,&(CG2gzi_vel),CG2s,&(CG2s_vel),
			     tau2,acc_TermCG2,acc_TermCG22,CG2vel_Term,TCG2,TCG2o);

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

    ABAb_integ_cort(qrotCG1,CG1qvel,CG1q,qaccCG1,CG1predict,CG1correct,dt,numclut);
    ABAb_integ_cort_Term(CG1predict_Term,CG1predict_Term2,
			 CG1correct_Term,CG1correct_Term2,
			 acc_TermCG1,acc_TermCG12,CG1vel_Term,
			 delta_TermCG1,dt);
    ABAbNH_update_cort_new(&CG1gzi,CG1gzi_vel,&CG1s,&CG1s_vel,
			   CG1predict_gzi,CG1correct_gzi,
			   CG1predict_s,CG1correct_s,dt);
    ABAb_update(CG1clt,CG1crd,qrotCG1,numclut,numatom);

    ABAb_integ_cort(qrotCG2,CG2qvel,CG2q,qaccCG2,CG2predict,CG2correct,dt,numclut);
    ABAb_integ_cort_Term(CG2predict_Term,CG2predict_Term2,
			 CG2correct_Term,CG2correct_Term2,
			 acc_TermCG2,acc_TermCG22,CG2vel_Term,
			 delta_TermCG2,dt);
    ABAbNH_update_cort_new(&CG2gzi,CG2gzi_vel,&CG2s,&CG2s_vel,
			   CG2predict_gzi,CG2correct_gzi,
			   CG2predict_s,CG2correct_s,dt);
    ABAb_update(CG2clt,CG2crd,qrotCG2,numclut,numatom);

    TACCM_integ_cort_Z(predict_Z,correct_Z,accZ,Z,velZ,numZ,dt,pi);
    TACCM_NH_update_cort_new(&gziZ,gzi_velZ,&sZ,&s_velZ,predict_gziZ,correct_gziZ,predict_sZ,correct_sZ,dt);

    ABAb_calcKineE_TermOn_new_simp(&KEFG,FGclt,FGcrd,FGqvel,FGvel_Term,numclut,numatom);
    ABAbNH_calcKE_new(FGgzi,FGs,FGs_vel,FGQ,KEFGo,&PEvFG,&KEvFG);

    ABAb_calcKineE_TermOn_new_simp(&KECG1,CG1clt,CG1crd,CG1qvel,CG1vel_Term,numclut,numatom);
    ABAbNH_calcKE_new(CG1gzi,CG1s,CG1s_vel,CG1Q,KECG1o,&PEvCG1,&KEvCG1);

    ABAb_calcKineE_TermOn_new_simp(&KECG2,CG2clt,CG2crd,CG2qvel,CG2vel_Term,numclut,numatom);
    ABAbNH_calcKE_new(CG2gzi,CG2s,CG2s_vel,CG2Q,KECG2o,&PEvCG2,&KEvCG2);

    TACCM_calcKineE_Z(&KEZ,massZ,velZ,numZ);
    TACCM_NH_calcKE_new(gziZ,sZ,s_velZ,/*Q_NH*/QZ,/*KEobj*/KEZo,&PEvZ,&KEvZ);
    TZ=KEZ/(numZ*k_B)*2.0;

    if (i%interval==0) {
      //      KEFG=KEFG/UNITT;     TFG=KEFG/(DOF*k_B)*2.0;
      //      PEvFG=PEvFG/UNITT;   KEvFG=KEvFG/UNITT;
      
      //      KECG1=KECG1/UNITT;     TCG1=KECG1/(DOF*k_B)*2.0;
      //      PEvCG1=PEvCG1/UNITT;  KEvCG1=KEvCG1/UNITT;

      PEFG=0.5*e.p_e_t+0.5*e.p_LJ_t+0.5*e.p_e_14_t+0.5*e.p_LJ_14_t+e.p_d_t/*+e.p_a_t+e.p_b_t*/;
      EtFG=PEFG+KEFG+PEvFG+KEvFG;
      fprintf(outputfileFG,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "
                  	   ,i+1,PEFG, KEFG, KEvFG,PEvFG,EtFG,*PEZFG,TFG);

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
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) COM[k]+=mass[j]*FGcrd[j*3+k]/summass;
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) FGcrd[j*3+k]-=COM[k];
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=FGcrd[j*3+k];
      myncL_put_crd_ene_MCD(nc_id_MCDFG,*l,crd_nc,e,0.0);

      summass=0.0; for (j=0;j<numatom;++j) summass+=mass[j];
      for (j=0;j<3;++j) COM[j]=0.0;
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) COM[k]+=mass[j]*CG1crd[j*3+k]/summass;
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) CG1crd[j*3+k]-=COM[k];
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=CG1crd[j*3+k];
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=CG1crd[j*3+k];
      myncL_put_crd_ene_MCD(nc_id_MCDCG1,*l,crd_nc,e,0.0);

      summass=0.0; for (j=0;j<numatom;++j) summass+=mass[j];
      for (j=0;j<3;++j) COM[j]=0.0;
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) COM[k]+=mass[j]*CG2crd[j*3+k]/summass;
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) CG2crd[j*3+k]-=COM[k];
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=CG2crd[j*3+k];
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=CG2crd[j*3+k];
      myncL_put_crd_ene_MCD(nc_id_MCDCG2,*l,crd_nc,e,0.0);

      ++(*l);

      ///////////////// TACCM //////////////////////
      //      KEZ=KEZ/UNITT;      TZ=KEZ/(numZ*k_B)*2.0;

      //      PEvZ=PEvZ/UNITT;      KEvZ=KEvZ/UNITT;

      EtZ=*PEZ+KEZ+PEvZ+KEvZ;
      fprintf(outputfileFG,"%d %e %e %e %e %e %e \n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);

      fprintf(outputfileCG1,"%d %e %e %e %e %e %e \n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);

      fprintf(outputfileCG2,"%d %e %e %e %e %e %e \n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);

      for (j=0;j<numZ;++j) fprintf(trjfileZ,"%e ",Z[j]);
      fprintf(trjfileZ,"\n");
      for (j=0;j<numZ;++j) fprintf(trjfilThetaFG,"%e ",thetaFG[j]);
      fprintf(trjfilThetaFG,"\n");      
      for (j=0;j<numZ;++j) fprintf(trjfilThetaCG1,"%e ",thetaCG1[j]);
      fprintf(trjfilThetaCG1,"\n");      
      for (j=0;j<numZ;++j) fprintf(trjfilThetaCG2,"%e ",thetaCG2[j]);
      fprintf(trjfilThetaCG2,"\n");      
      ///////////////// TACCM //////////////////////
    }
  }

  return *PEZ;
}
