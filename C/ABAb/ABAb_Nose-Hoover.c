
#include <stdio.h>
#include <math.h>

#include "ABAb.h"

void ABAbNH_update_pret(double *s,double *s_vel,double predict_s[6], double correct_s[6],double dt) {
  int i,j;

  for (i=0;i<6;++i) predict_s[i] = 0.0;
  for (i=0;i<6;++i) for (j=0;j<6;++j) predict_s[i] += Telar_Matrix[i][j]*correct_s[j];
  *s=predict_s[0];
  *s_vel=predict_s[1]/dt;
}
void ABAbNH_update_cort(double *gzi, double *s, double *s_vel, double s_acc, double predict_s[6], double correct_s[6], double dt) {
  int i,j,k;
  double ds;

  ds = 0.5*dt*dt*s_acc-predict_s[2];
  for (i=0;i<6;++i) correct_s[i] = predict_s[i]+GearsConstant[i]*ds;

  *s=correct_s[0];
  *s_vel=correct_s[1]/dt;
  *gzi=(*s_vel)/((*s)*(*s));
}

double ABAbNH_calcKE(double s,double s_vel,double Q,double KEobj,double *PEv,double *KEv){
  int i,j,k;

  *KEv = 0.5*Q*(s_vel/s)*(s_vel/s);
  *PEv = KEobj*log(s);
}

void ABAbNH_set(double s,double s_vel,double gzi,double predict_s[6],double correct_s[6],double tau, double *tau2, double *Q, double KEobj,double dt) {
  int i,j;
  double UNIT=418.4070;

  *tau2=tau*tau;
  //  *Q=(*tau2)*KEobj*2.0
  *Q=(*tau2)*KEobj/**2.0*/;

  for (i=0;i<6;++i){
    correct_s[i]=0.0;
    predict_s[i]=0.0;
  }

  correct_s[0]=s;
  correct_s[1]=dt*s_vel;

  gzi=s_vel/s;
}

void ABAbm_backpass_TermOn_NH(double *qacc,ABIb* abi,CLTb *clt,int numclut,double *qvel,double s,double s_vel,double *s_acc,double tau2,double *acc_Term,double *acc_Term2,double *vel_Term,double Temp,double TempB) {
  int i,j,k;
  double add[6];
  int nparent;
  double UNIT=418.4070;

  for(i=0;i<6;++i) clt[0].Spacc[i]=acc_Term2[i];
  for (i=1;i<numclut;++i) {
    nparent=clt[i].nNumClutOfParent-1;

    for(j=0;j<6;++j) {
      clt[i].Spacc[j] = 0.0;
      clt[i].preSpacc[j] = 0.0;
    }

    for (j=0;j<6;++j) 
      for (k=0;k<6;++k) 
	clt[i].preSpacc[j]+=clt[i].TM[k][j]*clt[nparent].Spacc[k];
    
    qacc[i] = 0.0;
    for(j=0;j<6;++j)
      qacc[i]-=abi[i].KG[j]*clt[i].preSpacc[j];
    qacc[i]+=abi[i].nyu;
    
    for(j=0;j<6;++j) clt[i].Spacc[j]+=clt[i].preSpacc[j]+clt[i].Coacc[j];
    clt[i].Spacc[2]+=qacc[i];
  }

  add[0]=0.0;
  add[1]=0.0;
  add[2]=0.0;
  add[3]=(-vel_Term[2]*vel_Term[3+1]+vel_Term[1]*vel_Term[3+2]);
  add[4]=( vel_Term[2]*vel_Term[3+0]-vel_Term[0]*vel_Term[3+2]);
  add[5]=(-vel_Term[1]*vel_Term[3+0]+vel_Term[0]*vel_Term[3+1]);

  for(i=0;i<6;++i) { 
    acc_Term2[i] -= s_vel/(s)*vel_Term[i];
    acc_Term[i]=acc_Term2[i]-add[i];
  }

  for(i=0;i<numclut;++i) qacc[i] -= s_vel/(s)*qvel[i];
  //  *s_acc = 1.0/(tau2)*s*(2.0*Temp/TempB-1.0)+s_vel*s_vel/s;
  *s_acc = 1.0/(tau2)*s*(Temp/TempB-1.0)+s_vel*s_vel/s;
}

void ABAbm_backpass_NH(double *qacc,ABIb* abi,CLTb *clt,int numclut,double *qvel,double s,double s_vel,double *s_acc,double tau2,double Temp,double TempB) {
  int i,j,k;
  int nparent;
  double UNIT=418.4070;

  for(i=0;i<6;++i) clt[0].Spacc[i]=0.0;
  for (i=1;i<numclut;++i) {
    nparent=clt[i].nNumClutOfParent-1;

    for(j=0;j<6;++j) {
      clt[i].Spacc[j] = 0.0;
      clt[i].preSpacc[j] = 0.0;
    }

    for (j=0;j<6;++j) 
      for (k=0;k<6;++k) 
	clt[i].preSpacc[j]+=clt[i].TM[k][j]*clt[nparent].Spacc[k];
    
    qacc[i] = 0.0;
    for(j=0;j<6;++j)
      qacc[i]-=abi[i].KG[j]*clt[i].preSpacc[j];
    qacc[i]+=abi[i].nyu;
    
    for(j=0;j<6;++j) clt[i].Spacc[j]+=clt[i].preSpacc[j]+clt[i].Coacc[j];
    clt[i].Spacc[2]+=qacc[i];
  }

  for(i=0;i<numclut;++i) qacc[i] -= s_vel/(s)*qvel[i];
  //  *s_acc = 1.0/(tau2)*s*(2.0*Temp/TempB-1.0)+s_vel*s_vel/s;
  *s_acc = 1.0/(tau2)*s*(Temp/TempB-1.0)+s_vel*s_vel/s;
}

void ABAbm_corBF_NH_TERM(double corABIb[6][6],double corBF[6],double *acc_Term,double *acc_Term2,double *vel_Term,double s,double s_vel){
  int i,j,k;
  double *temp,*inv,add[6];

  temp=(double *)gcemalloc(sizeof(double)*6*6);
  inv=(double *)gcemalloc(sizeof(double)*6*6);

  for (i=0;i<6;++i)
    for (j=0;j<6;++j)
      temp[i*6+j]=corABIb[i][j];
  invm2(temp,inv,6);
    for(i=0;i<6;++i) {
      acc_Term[i]=0.0;
      acc_Term2[i]=0.0;
    }
    add[0]=0.0;
    add[1]=0.0;
    add[2]=0.0;
    add[3]=(-vel_Term[2]*vel_Term[3+1]+vel_Term[1]*vel_Term[3+2]);
    add[4]=( vel_Term[2]*vel_Term[3+0]-vel_Term[0]*vel_Term[3+2]);
    add[5]=(-vel_Term[1]*vel_Term[3+0]+vel_Term[0]*vel_Term[3+1]);

    for(i=0;i<6;++i)
      for(j=0;j<6;++j)
	acc_Term2[i]-=inv[i*6+j]*corBF[j];
}

void ABAbm_mainpass_TermOn_NH(ABIb* abi,CLTb *clt,double *Q,int numclut,double *acc_Term,double *acc_Term2,double *vel_Term,double s,double s_vel) {
  int i,j,k,nb;
  int nparent;
  double TMI[4][6][6],preABIbI[4][6][6],preBFI[4][6];

  for (i=numclut-1;i>=0;--i) {
    if (clt[i].terminal==TERM) {
       for (nb=0;nb<4;++nb) {
	for (j=0;j<6;++j) {
	  preBFI[nb][j]=0.0;
	  for (k=0;k<6;++k) {
	    TMI[nb][j][k]=0.0;
	    preABIbI[nb][j][k]=0.0;
	  }
	}
      }
      ABAbm_cCABIb(abi[i].corABIb,preABIbI,clt[i].IM,TMI,clt[i].num_branch);
    }
    else {
      for (nb=0;nb<clt[i].num_branch;++nb) {
	for (j=0;j<6;++j) {
	  for (k=0;k<6;++k) {
	    TMI[nb][j][k]=clt[clt[i].nNumClutOfChild[nb]-1].TM[j][k];
	    preABIbI[nb][j][k]=abi[clt[i].nNumClutOfChild[nb]-1].preABIb[j][k];
	  }
	}
      }
      ABAbm_cCABIb(abi[i].corABIb,preABIbI,clt[i].IM,TMI,clt[i].num_branch);
    }

    abi[i].D=ABAbm_cD(abi[i].corABIb);
    ABAbm_cKG(abi[i].KG,abi[i].corABIb);
    ABAbm_cPABIb(abi[i].preABIb,abi[i].corABIb);

    if (clt[i].terminal==TERM) {
      for (nb=0;nb<4;++nb) {
	for (j=0;j<6;++j) {
	  preBFI[nb][j]=0.0;
	  for (k=0;k<6;++k) {
	    TMI[nb][j][k]=0.0;
	  }
	}
      }
      ABAbm_corBF(abi[i].corBF,preBFI,abi[i].corABIb,clt[i].Coacc,clt[i].Cofrc,clt[i].Spfrc,TMI,clt[i].num_branch);
    }
    else {
      for (nb=0;nb<clt[i].num_branch;++nb) {
	for (j=0;j<6;++j) {
	  preBFI[nb][j]=abi[clt[i].nNumClutOfChild[nb]-1].preBF[j];
	  for (k=0;k<6;++k) {
	    TMI[nb][j][k]=clt[clt[i].nNumClutOfChild[nb]-1].TM[j][k];
	  }
	}
      }
      ABAbm_corBF(abi[i].corBF,preBFI,abi[i].corABIb,clt[i].Coacc,clt[i].Cofrc,clt[i].Spfrc,TMI,clt[i].num_branch);
    }
    abi[i].eata=ABAbm_eata(Q[i],abi[i].corBF);
    abi[i].nyu=ABAbm_nyu(abi[i].eata,abi[i].D);
    ABAbm_preBF(abi[i].preBF,abi[i].corBF,abi[i].KG,abi[i].eata);
  }
  ABAbm_corBF_NH_TERM(abi[0].corABIb,abi[0].corBF,acc_Term,acc_Term2,vel_Term,s,s_vel);
}


