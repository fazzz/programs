
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "ABA.h" // 2014-06-18
#include "ABAb.h"  // 2014-06-18
#include "EF.h"
#include "LA.h"

void ABAp_prepass(CLT *clt,double *qvel,int numclt,int numatom,double *crd) {
  int i,j;
  int nNumClutOfParent;

  ABAs_trans_Matrix(clt,numclt,numatom,crd);

  for (i=0;i<numclt;++i) {
    nNumClutOfParent = clt[i].nNumClutOfParent-1;

    if (i == 0) {
      for(j=0;j<6;++j) {
	clt[i].Spvel[j] = 0.0;
	clt[i].Cofrc[j] = 0.0;
	clt[i].Coacc[j] = 0.0;
      }
    }
    else {
      ABAp_Spvelo(clt[i].Spvel,qvel[i],clt[i].TM,clt[nNumClutOfParent].Spvel);
    }
  }

  for (i=0;i<numclt;++i) {
    nNumClutOfParent = clt[i].nNumClutOfParent-1;
    ABAp_Cofr(clt[i].Cofrc,clt[i].Spvel,clt[i].IM,clt[i].qCOM,clt[i].sum_mass);
    ABA_Coacc(clt[i].Coacc,clt[i].Spvel,clt[nNumClutOfParent].Spvel,qvel[i],
	      clt[i].origin_atom_a,clt[nNumClutOfParent].origin_atom_a,nNumClutOfParent,
	      crd,clt[nNumClutOfParent].trans_A_to_CN,clt[i].TM);
  }
}

void ABAp_prepass_TermOn(CLT *clt,double *qvel,int numclt,int numatom,double *crd,double *vel_Term) {
  int i,j;
  int nNumClutOfParent;

  ABAs_trans_Matrix(clt,numclt,numatom,crd);

  for (i=0;i<numclt;++i) {
    nNumClutOfParent = clt[i].nNumClutOfParent-1;

    if (i == 0) {
      for(j=0;j<6;++j) {
	clt[i].Spvel[j] = vel_Term[j];
	clt[i].Cofrc[j] = 0.0;
	clt[i].Coacc[j] = 0.0;
      }
    }
    else {
      ABAp_Spvelo(clt[i].Spvel,qvel[i],clt[i].TM,clt[nNumClutOfParent].Spvel);
    }
  }

  for (i=0;i<numclt;++i) {
    nNumClutOfParent = clt[i].nNumClutOfParent-1;
    ABAp_Cofr(clt[i].Cofrc,clt[i].Spvel,clt[i].IM,clt[i].qCOM,clt[i].sum_mass);
    ABA_Coacc(clt[i].Coacc,clt[i].Spvel,clt[nNumClutOfParent].Spvel,qvel[i],
	      clt[i].origin_atom_a,clt[nNumClutOfParent].origin_atom_a,nNumClutOfParent,
	      crd,clt[nNumClutOfParent].trans_A_to_CN,clt[i].TM);
  }
}

void ABAp_Cofr(double Cofrc[6],double Spvel[6],double IM[6][6],double qCOM[3],double Sum_Mass) {
  int i,j,k;
  int alpha, alpha2,alpha3;

  double omega[3];
  double Inertia[3][3];
  
  double omegaproduct[3][3];
  double MasproductS[3];
  double matOmegaOmega[3][3];
  double InertiaProductOmega[3];
  double MasproductSProductOmega[3];
  
  double omegaqCOMvelo[3];
  double omegaqCOM[3][3];
  double qCOMproduct[3][3];
  

  for (i=0;i<6;++i) Cofrc[i] = 0.0;

  for (alpha=0;alpha<3;++alpha) {
    InertiaProductOmega[alpha] = 0.0;
    MasproductSProductOmega[alpha] = 0.0;
  }

  omegaproduct[0][0]= 0.0;
  omegaproduct[0][1]=-Spvel[2];
  omegaproduct[0][2]= Spvel[1];
  omegaproduct[1][0]= Spvel[2];
  omegaproduct[1][1]= 0.0;
  omegaproduct[1][2]=-Spvel[0];
  omegaproduct[2][0]=-Spvel[1];
  omegaproduct[2][1]= Spvel[0];
  omegaproduct[2][2]= 0.0;

  for (alpha=0;alpha<3;++alpha)
    for (alpha2=0;alpha2<3;++alpha2)
      InertiaProductOmega[alpha]+=IM[alpha][alpha2]*Spvel[alpha2];

  for (alpha=0;alpha<3;++alpha)
    for (alpha2=0;alpha2<3;++alpha2)
      Cofrc[alpha]+=omegaproduct[alpha][alpha2]*InertiaProductOmega[alpha2];

  MasproductS[0] = qCOM[0]*Sum_Mass;
  MasproductS[1] = qCOM[1]*Sum_Mass;
  MasproductS[2] = qCOM[2]*Sum_Mass;

  for (alpha=0;alpha<3;++alpha)
    for (alpha2=0;alpha2<3;++alpha2)
      MasproductSProductOmega[alpha]
	+= omegaproduct[alpha][alpha2]*MasproductS[alpha2];

  for (alpha=0;alpha<3;++alpha)
    for (alpha2=0;alpha2<3;++alpha2)
      matOmegaOmega[alpha][alpha2] = 0.0;

  for (alpha=0;alpha<3;++alpha)
    for (alpha2=0;alpha2<3;++alpha2)
      for (alpha3=0;alpha3<3;++alpha3)
	matOmegaOmega[alpha][alpha2] += omegaproduct[alpha][alpha3]*omegaproduct[alpha3][alpha2];

  for (alpha=0;alpha<3;++alpha)
    for (alpha2=0;alpha2<3;++alpha2)
      Cofrc[alpha+3] += matOmegaOmega[alpha][alpha2]*MasproductS[alpha2];

}

void ABA_Coacc(double Coacc[6],double Spvel[6],double Spvel_P[6],double qvel,
	       int origin_atom_a,int origin_atom_a_Parent,int nNumminousone,
	       double *crd,double trans_A_to_CN_P[3][3],double TMat[6][6]) {
  int i,j;
  int alpha;
  int alpha2;
  int alpha3;
  int nNumAtomOfClut;
  int nNumAtomOfClutminousone;

  double RotatnNumtonNumMiOn[3][3];

  double Omega[3];
  double ez[3];

  double vector[3];
  double vector2[3];

  double Coriolis_acc_dummy_up[3];
  double Coriolis_acc_dummy_up2[3];
  double Coriolis_acc_dummy[3];
  double matOmegaOmega[3][3];
  double omegaproduct[3][3];
  double omegaproduct2[3][3];
  double Coriolis_acc_dummy_dummy[3];

  for (i=0;i<6;++i) Coacc[i] = 0.0;

  for (alpha=0;alpha<3;++alpha) {
    Omega[alpha] = 0.0;
    vector[alpha] = 0.0;
    vector2[alpha] = 0.0;
    Coriolis_acc_dummy[alpha] = 0.0;
    Coriolis_acc_dummy_up[alpha] = 0.0;
    Coriolis_acc_dummy_up2[alpha] = 0.0;
  }

  for (alpha=0;alpha<3;++alpha) Omega[alpha]= Spvel[alpha];

  ez[0] = 0.0;
  ez[1] = 0.0;
  ez[2] = 1.0;
  
  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      Coriolis_acc_dummy_up[i] = Omega[i];

  omegaproduct2[0][0]= 0.0;
  omegaproduct2[0][1]=-Coriolis_acc_dummy_up[2];
  omegaproduct2[0][2]= Coriolis_acc_dummy_up[1];
  omegaproduct2[1][0]= Coriolis_acc_dummy_up[2];
  omegaproduct2[1][1]= 0.0;
  omegaproduct2[1][2]=-Coriolis_acc_dummy_up[0];
  omegaproduct2[2][0]=-Coriolis_acc_dummy_up[1];
  omegaproduct2[2][1]= Coriolis_acc_dummy_up[0];
  omegaproduct2[2][2]= 0.0;

  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      Coriolis_acc_dummy_up2[i]+= omegaproduct2[i][j]*ez[j];

  for (i=0;i<3;++i) Coacc[i] = Coriolis_acc_dummy_up2[i]*qvel;

  nNumAtomOfClut = origin_atom_a-1;
  if(nNumminousone==-1) 
    for (alpha=0;alpha<3;++alpha)
      vector[alpha] = crd[nNumAtomOfClut*3+alpha];
  else {
    nNumAtomOfClutminousone = origin_atom_a_Parent-1;
    for (alpha=0;alpha<3;++alpha)
      vector[alpha] = crd[nNumAtomOfClut*3+alpha]-crd[nNumAtomOfClutminousone*3+alpha];
  }
  
  for (alpha=0;alpha<3;++alpha) vector2[alpha] = 0.0;
  
  for (alpha=0;alpha<3;++alpha)
    for (alpha2=0;alpha2<3;++alpha2)
      vector2[alpha] += trans_A_to_CN_P[alpha][alpha2]*vector[alpha2];

  omegaproduct[0][0]= 0.0;
  omegaproduct[0][1]=-Spvel_P[2];
  omegaproduct[0][2]= Spvel_P[1];
  omegaproduct[1][0]= Spvel_P[2];
  omegaproduct[1][1]= 0.0;
  omegaproduct[1][2]=-Spvel_P[0];
  omegaproduct[2][0]=-Spvel_P[1];
  omegaproduct[2][1]= Spvel_P[0];
  omegaproduct[2][2]= 0.0;

  for (alpha=0;alpha<3;++alpha) {
    Coriolis_acc_dummy[alpha] = 0.0;
    Coriolis_acc_dummy_dummy[alpha] = 0.0;
  }

  for (alpha=0;alpha<3;++alpha)
    for (alpha2=0;alpha2<3;++alpha2)
      Coriolis_acc_dummy_dummy[alpha] += omegaproduct[alpha][alpha2]*vector2[alpha2];

  for (alpha=0;alpha<3;++alpha)
    for (alpha2=0;alpha2<3;++alpha2)
      Coriolis_acc_dummy[alpha] += omegaproduct[alpha][alpha2]*Coriolis_acc_dummy_dummy[alpha2];

  for (alpha=0;alpha<3;++alpha)
    for (alpha2=0;alpha2<3;++alpha2)
      Coacc[alpha+3]+= TMat[alpha2][alpha]*Coriolis_acc_dummy[alpha2];
}

void ABAp_Spvelo(double Spvel[6],double qvel,double TMat[6][6],double SpvelPart[6]) {
  int i,j;

  for(i=0;i<6;++i) Spvel[i] = 0.0;

  for(i=0;i<6;++i)
    for(j=0;j<6;++j)
      Spvel[i]+=TMat[j][i]*SpvelPart[j];

  Spvel[2]+=qvel;

}

