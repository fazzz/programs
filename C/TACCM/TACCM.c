
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MD_NHC_MP1996.h"
#include "FFL.h"
//#include "UMBP.h"
#include "PTL.h"

#include "TOPO.h"
#include "LA.h"
#include "mymath.h"

#include "ABAb.h"
#include "TACCM.h"

#include "EF.h"

#include "RAND.h"
#include "BOXMULL.h"

#define UNITT 418.4070

#define ON  1
#define OFF 0

//#define nys 3

/******************************************************************************************************************************************************************/
/* double TACCM_calc_eff_f_dihed_type_AAFF_Amber(double *crd,int numatom,											  */
/* 					      double *z,  int numz,												  */
/* 					      int *pairp,double *fcp,												  */
/* 					      double **f_crd, double *f_z,											  */
/* 					      struct potential *e, struct force *f){										  */
/*   int i,j,k;																			  */
/*   ndouble *p_UMB,p_UMB_t,**f_UMB,*theta;															  */
/*   dnouble p_t;																		  */
/* 																				  */
/*   ffL_calcffandforce(crd,numatom,p,f);															  */
/*   theta=UMB_calc_dihetype_ff(crd,numatom,pairp,nump,fcp,z,p_UMB,f_UMB);											  */
/* 																				  */
/*   for(i=0;i<numatom;++i)																	  */
/*     for(j=0;j<3;++j)																		  */
/*       f_crd[i][j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]+(*f).f_e_14[i*3+j]+(*f).f_LJ_14[i*3+j]+(*f).f_e[i*3+j]+(*f).f_LJ[i*3+j]+f_UMB[i][j];	  */
/* 																				  */
/*   for(i=0;i<nump;++i) {																	  */
/*     f_z[i]=fcp*(theta[i]-z[i]);																  */
/*     p_UMB_t+=p_UMB[i];																	  */
/*   }																				  */
/* 																				  */
/*   p_t=0.5*e.p_e_t+0.5*e.p_LJ_t+0.5*e.p_e_14_t+0.5*e.p_LJ_14_t+e.p_d_t+e.p_a_t+e.p_b_t+pUMB_t;								  */
/* 																				  */
/*   return p_t;																		  */
/* }																				  */
/******************************************************************************************************************************************************************/

double TACCM_calc_eff_FF(double *theta, double *Z,  int numZ,double Kapa, double *Q, int **pairs,double pi){
  int i,j,k;
  double PE=0.0;
  double delta;

  for (i=0;i<numZ;++i) {
    //    delta=Z[i]-theta[i];
    //    if (delta < 0) delta=-1.0*delta;
    //    if (delta>pi) delta=2.0*pi-delta;

    if ((delta=Z[i]-theta[i])>pi) delta-=2.0*pi;
    else if ((delta=Z[i]-theta[i])<-1.0*pi) delta+=2.0*pi;

    PE+=0.5*Kapa*delta*delta;
    Q[pairs[i][4]]=-Kapa*delta*UNITT;
  }

  return PE;
}

double TACCM_calc_eff_FF_Z(double *Z,int numZ,double *theta,double KZ,double *f, double pi){
  int i,j;
  double PE=0.0;
  double delta;

  for (i=0;i<numZ;++i) {
    //    delta=Z[i]-theta[i];
    //    if (delta < 0) delta=-1.0*delta;
    //    if (delta>pi) delta=2.0*pi-delta;
    if ((delta=Z[i]-theta[i])>pi) delta-=2.0*pi;
    else if ((delta=Z[i]-theta[i])<-1.0*pi) delta+=2.0*pi;
    PE+=0.5*KZ*delta*delta;
    f[i]=-KZ*delta*UNITT;
  }

  return PE;
}

double TACCM_calc_eff_FF_Z_2(double *Z,int numZ,double *theta,double KZ,double *f, double *PE,double pi){
  int i,j;
  //  double PE=0.0;
  double delta;

  *PE=0.0;
  for (i=0;i<numZ;++i) {
    //    delta=Z[i]-theta[i];
    //    if (delta < 0) delta=-1.0*delta;
    //    if (delta>pi) delta=2.0*pi-delta;
    if ((delta=Z[i]-theta[i])>pi) delta-=2.0*pi;
    else if ((delta=Z[i]-theta[i])<-1.0*pi) delta+=2.0*pi;
    *PE+=0.5*KZ*delta*delta;
    f[i]=-KZ*delta*UNITT;
  }

  return *PE;
}

double TACCM_CTheta(double *crd,int numatom,double *theta, int numdihe, int **pairs, double pi){
  int i,j,k,l;
  int ii,jj,kk,ll;

  double m[3],n[3],m_n[3],n_n[3],lm,ln;
  double vij[3],vkj[3],vkl[3];
  double lkj;
  double vijvkj,vklvkj;

  double atom[4][3];
  double dihed;

  for (i=0;i<numdihe;++i) {
    ii=pairs[i][0]-1;
    jj=pairs[i][1]-1;
    kk=pairs[i][2]-1;
    ll=pairs[i][3]-1;

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

  return 0.0;
}

double TACCM_MD_Generate_inivelo(double *velZ,double massZ,int numZ,double KbT) {
  int i,j;
  double KE=0.0;
  //  double UNITT=418.4070;
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  for (i=0;i<numZ;++i) velZ[i]=Box_Muller(i,0.0,KbT/massZ*0.5);

  for (i=0;i<numZ;++i) KE+=0.5*massZ*velZ[i]*velZ[i];

  return KE;
}

void TACCM_integ_pret_Z(double **predict_Z,double **correct_Z,double *Z,double *vel_Z,int numZ,double dt,double pi){
  int i,j,k;

  for (i=0;i<numZ;++i) 
    for (j=0;j<6;++j) 
      predict_Z[i][j] = 0.0;

  for (i=0;i<numZ;++i) 
    for (j=0;j<6;++j) 
      for (k=0;k<6;++k) 
  	predict_Z[i][j] += Telar_MatrixZ[j][k]*correct_Z[i][k];

  // debug //
  for (i=0;i<numZ;++i) {
    while (predict_Z[i][0]<-1.0*pi) predict_Z[i][0]+=2.0*pi;
    while (predict_Z[i][0]>pi) predict_Z[i][0]-=2.0*pi;
  }
  // debug //

  for (i=0;i<numZ;++i) {
    Z[i]=predict_Z[i][0];
    vel_Z[i]=predict_Z[i][1]/dt;
  }
  
}

void TACCM_integ_cort_Z(double **predict_Z,double **correct_Z,double *acc_Z,double *Z,double *vel_Z,int numZ,double dt,double pi) {
  int i,j,k;
  double *delta_acc_Z;

  delta_acc_Z=(double *)gcemalloc(sizeof(double)*numZ);
  
  for (i=0;i<numZ;++i) delta_acc_Z[i] = 0.5*dt*dt*acc_Z[i]-predict_Z[i][2];

  for (i=0;i<numZ;++i) for (j=0;j<6;++j) correct_Z[i][j] = 0.0;

  for (i=0;i<numZ;++i) 
    for (j=0;j<6;++j)   
      correct_Z[i][j] = predict_Z[i][j]+GearsConstantZ[j]*delta_acc_Z[i];

  // debug //
  for (i=0;i<numZ;++i) {
    while (correct_Z[i][0]<-1.0*pi) correct_Z[i][0]+=2.0*pi;
    while (correct_Z[i][0]>pi) correct_Z[i][0]-=2.0*pi;
  }
  // debug //
  
  for (i=0;i<numZ;++i) {
    Z[i]=correct_Z[i][0];
    vel_Z[i]=correct_Z[i][1]/dt;
  }
}

double TACCM_calcKineE_Z(double *KE,double massZ,double *vel_Z,int numZ) {
  int i,j,k;

  *KE=0.0;

  for (i=0;i<numZ;++i) {
    *KE+=massZ*vel_Z[i]*vel_Z[i];
  }
  *KE=0.5*(*KE)/(4.18407*100.0)/*/UNIT*/;

}

double TACCM_solver_NH_Z(double *accZ,double *velZ,double massZ,double *frcZ,int numZ,double gzi,double *gzi_vel,double tau2,double Temp,double TempB){
  int i,j,k;

  for (i=0;i<numZ;++i) {
    accZ[i]=frcZ[i]/massZ-gzi*velZ[i];
  }

  //  for(i=0;i<numZ;++i) { 
  //    accZ[i] -= gzi*velZ[i];
  //  }

  *gzi_vel = 1.0/(tau2)*(Temp/TempB-1.0);

}

void TACCM_NH_update_pret_new(double *gzi,double *gzi_vel,double predict_gzi[5], double correct_gzi[5],
			      double *s,double *s_vel,double predict_s[5], double correct_s[5],double dt) {
  int i,j;

  for (i=0;i<5;++i) {
    predict_gzi[i] = 0.0;
    predict_s[i] = 0.0;
  }

  for (i=0;i<5;++i) {
    for (j=0;j<5;++j) { 
      predict_gzi[i] += Telar_MatrixZ[i][j]*correct_gzi[j];
      predict_s[i]   += Telar_MatrixZ[i][j]*correct_s[j];
    }
  }

  *gzi=predict_gzi[0];
  *s=predict_s[0];
  *s_vel=predict_s[1]/dt;
}

void TACCM_NH_update_cort_new(double *gzi, double gzi_vel, double *s, double *s_vel,
			      double predict_gzi[5], double correct_gzi[5],
			      double predict_s[5], double correct_s[5], double dt) {
  int i,j,k;
  double d,d2;

  d = dt*gzi_vel-predict_gzi[1];
  for (i=0;i<5;++i) correct_gzi[i] = predict_gzi[i]+GearsConstantZ5[i]*d;
  *gzi=correct_gzi[0];

  *s_vel   = *gzi*(*s)/**(*s)*/;
  d2 = dt*(*s_vel)-predict_s[1];
  for (i=0;i<5;++i) correct_s[i] = predict_s[i]+GearsConstantZ5[i]*d2;
  *s=correct_s[0];
  *s_vel=correct_s[1]/dt;

}

double TACCM_NH_calcKE_new(double gzi,double s,double s_vel,double Q,double KEobj,double *PEv,double *KEv){
  int i,j,k;

  *KEv = 0.5*Q*(s_vel/s)*(s_vel/s);
  *PEv = 2.0*KEobj*log(s);
}

void TACCM_NH_set_new(double s,double s_vel,double gzi,double predict_gzi[5],double correct_gzi[5],
		      double predict_s[5],double correct_s[5],double tau, double *tau2, 
		      double *Q, double KEobj,double dt) {
  int i,j;
  double pi;

  pi=acos(-1.0);

  tau=tau/2.0/pi;
  *tau2=tau*tau;

  *Q=(*tau2)*KEobj*2.0;   

  for (i=0;i<5;++i){
    correct_gzi[i]=0.0;
    predict_gzi[i]=0.0;
    correct_s[i]=0.0;
    predict_s[i]=0.0;
  }

  GearsConstantZ5[0]=250.0/720.0;
  GearsConstantZ5[1]=1.0;
  GearsConstantZ5[2]=11.0/12.0;
  GearsConstantZ5[3]=1.0/3.0;
  GearsConstantZ5[4]=1.0/24.0;

  correct_s[0]=s;
  correct_s[1]=dt*s_vel;

  gzi=s_vel/s;

  correct_gzi[0]=gzi;

  GearsConstantZ[0] = 3.0/16.0;
  GearsConstantZ[1] = 251.0/360.0;
  GearsConstantZ[2] = 1.0;
  GearsConstantZ[3] = 11.0/18.0;
  GearsConstantZ[4] = 1.0/6.0;
  GearsConstantZ[5] = 1.0/60.0;

  for (i=0;i<6;++i){
    for (j=0;j<6;++j){
      if (i != j) Telar_MatrixZ[i][j] = 0.0;
      else Telar_MatrixZ[i][i] = 1.0;
    }
  }
	
  Telar_MatrixZ[0][1] = 1.0;
  Telar_MatrixZ[0][2] = 1.0;
  Telar_MatrixZ[0][3] = 1.0;
  Telar_MatrixZ[0][4] = 1.0;
  Telar_MatrixZ[0][5] = 1.0;
  Telar_MatrixZ[1][2] = 2.0;
  Telar_MatrixZ[1][3] = 3.0;
  Telar_MatrixZ[1][4] = 4.0;
  Telar_MatrixZ[1][5] = 5.0;
  Telar_MatrixZ[2][3] = 3.0;
  Telar_MatrixZ[2][4] = 6.0;
  Telar_MatrixZ[2][5] = 10.0;
  Telar_MatrixZ[3][4] = 4.0;
  Telar_MatrixZ[3][5] = 10.0;
  Telar_MatrixZ[4][5] = 5.0;
}
