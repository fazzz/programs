
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "DCA.h"
#include "EF.h"
#include "LA.h"

void DCAp_prepass(CLT *clt,double *qvel,int numclt,int numatom) {
  int i;
  int nNumClutOfParent;

  for (i=0;i<numclt;++i)

  nNumClutOfParent = clt[nNumClut].nNumClutOfParent-1;

  if (nNumClut == 0) for(i=0;i<6;++i) clt[nNumClut].Coriolis_b[i] = 0.0;
  else sub_set_coriolis_force(clt,nNumClut,nNumClutOfParent);

}

void sub_set_coriolis_force(CLT *clt,int nNumClut, int nNumClutminousone) {
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
  

  for (i=0;i<6;++i) clt[nNumClut].Coriolis_b[i] = 0.0;

  for (alpha=0;alpha<3;++alpha) {
    InertiaProductOmega[alpha] = 0.0;
    MasproductSProductOmega[alpha] = 0.0;
  }

  omegaproduct[0][0]= 0.0;
  omegaproduct[0][1]=-clt[nNumClut].sp_velo[2];
  omegaproduct[0][2]= clt[nNumClut].sp_velo[1];
  omegaproduct[1][0]= clt[nNumClut].sp_velo[2];
  omegaproduct[1][1]= 0.0;
  omegaproduct[1][2]=-clt[nNumClut].sp_velo[0];
  omegaproduct[2][0]=-clt[nNumClut].sp_velo[1];
  omegaproduct[2][1]= clt[nNumClut].sp_velo[0];
  omegaproduct[2][2]= 0.0;

  for (alpha=0;alpha<3;++alpha)
    for (alpha2=0;alpha2<3;++alpha2)
      InertiaProductOmega[alpha]
	+= clt[nNumClut].InertiaMatrix[alpha][alpha2]*clt[nNumClut].sp_velo[alpha2];

  for (alpha=0;alpha<3;++alpha)
    for (alpha2=0;alpha2<3;++alpha2)
      clt[nNumClut].Coriolis_b[alpha]+=omegaproduct[alpha][alpha2]*InertiaProductOmega[alpha2];

  MasproductS[0] = clt[nNumClut].qCOM[0]*Sum_Mass(nNumClut);
  MasproductS[1] = clt[nNumClut].qCOM[1]*Sum_Mass(nNumClut);
  MasproductS[2] = clt[nNumClut].qCOM[2]*Sum_Mass(nNumClut);

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
      clt[nNumClut].Coriolis_b[alpha+3]
	+= matOmegaOmega[alpha][alpha2]*MasproductS[alpha2];

}

void set_coriolis_acc(CLT *clt,double **crd,int nNumClut) {
  int i;
  int nNumClutOfParent;

  nNumClutOfParent = clt[nNumClut].nNumClutOfParent-1;

  if (nNumClut == 0) for(i=0;i<6;++i) clt[nNumClut].Coriolis_acc[i] = 0.0;
  else sub_set_coriolis_acc(clt,crd,nNumClut, nNumClutOfParent);
  
}

void sub_set_coriolis_acc(CLT *clt,double **crd,int nNumClut,int nNumClutminusone) {
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

  for (i=0;i<6;++i) clt[nNumClut].Coriolis_acc[i] = 0.0;

  for (alpha=0;alpha<3;++alpha) {
    Omega[alpha] = 0.0;
    vector[alpha] = 0.0;
    vector2[alpha] = 0.0;
    Coriolis_acc_dummy[alpha] = 0.0;
    Coriolis_acc_dummy_up[alpha] = 0.0;
    Coriolis_acc_dummy_up2[alpha] = 0.0;
  }

  for (alpha=0;alpha<3;++alpha) Omega[alpha]= clt[nNumClut].sp_velo[alpha];

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

  for (i=0;i<3;++i) clt[nNumClut].Coriolis_acc[i] = Coriolis_acc_dummy_up2[i]*clt[nNumClut].ddihedang[0];

  nNumAtomOfClut = clt[nNumClut].origin_atom_a-1;
  if (nNumClutminusone==-1) vector[alpha] = crd[nNumAtomOfClut][alpha];
  else {
    nNumAtomOfClutminousone = clt[nNumClutminusone].origin_atom_a-1;
    for (alpha=0;alpha<3;++alpha)
      vector[alpha] = crd[nNumAtomOfClut][alpha]-crd[nNumAtomOfClutminousone][alpha];
  }
  
  for (alpha=0;alpha<3;++alpha) vector2[alpha] = 0.0;
  
  for (alpha=0;alpha<3;++alpha)
    for (alpha2=0;alpha2<3;++alpha2)
      vector2[alpha] += clt[nNumClutminusone].trans_A_to_CN[0][alpha][alpha2]*vector[alpha2];

  omegaproduct[0][0]= 0.0;
  omegaproduct[0][1]=-clt[nNumClutminusone].sp_velo[2];
  omegaproduct[0][2]= clt[nNumClutminusone].sp_velo[1];
  omegaproduct[1][0]= clt[nNumClutminusone].sp_velo[2];
  omegaproduct[1][1]= 0.0;
  omegaproduct[1][2]=-clt[nNumClutminusone].sp_velo[0];
  omegaproduct[2][0]=-clt[nNumClutminusone].sp_velo[1];
  omegaproduct[2][1]= clt[nNumClutminusone].sp_velo[0];
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
      clt[nNumClut].Coriolis_acc[alpha+3]+= clt[nNumClut].TransMatrix[0][alpha2][alpha]*Coriolis_acc_dummy[alpha2];
}

void set_sp_velo(CLT *clt,int nNumClut, int nNumClutOrigBranch){
  int i;
  int nNumClutOfParent;

  nNumClutOfParent = clt[nNumClut].nNumClutOfParent-1;

  if (nNumClut == 0)for(i=0;i<6;++i)clt[0].sp_velo[i] = 0.0;
  else sub_set_sp_velo(clt,nNumClut, nNumClutOfParent);

}

void sub_set_sp_velo(CLT *clt,int nNumClut, int nNumClutminusone) {
  int i,j,k;
  int alpha;
  int alpha2;
  int nNumFreedom;

  for(i=0;i<6;++i) clt[nNumClut].sp_velo[i] = 0.0;

  for(i=0;i<6;++i)
    for(j=0;j<6;++j)
      clt[nNumClut].sp_velo[i]+=clt[nNumClut].TransMatrix[0][j][i]*clt[nNumClutminusone].sp_velo[j];

  nNumFreedom = 2;

  clt[nNumClut].sp_velo[nNumFreedom]+= clt[nNumClut].ddihedang[0];

}

void set_trans_Matrix(CLT *clt,double **crd,int nNumClt,int nNumClutOrigBranch){
  int alpha;
  int alpha2;
  
  int nNumClutParent;

  nNumClutParent = clt[nNumClt].nNumClutOfParent-1;

  if (nNumClt == 0)
    for (alpha=0;alpha<6;++alpha)
      for (alpha2=0;alpha2<6;++alpha2)
	clt[nNumClt].TransMatrix[0][alpha][alpha2]= 0.0;
  else sub_set_trans_Matrix(clt,crd,nNumClt, nNumClutParent);

}

void sub_set_trans_Matrix(CLT *clt,double *crd,int nNumClt,int nNumCltminousone){
  int alpha,alpha2,alpha3,i,j,k;
  int nNumAtomOfClut,nNumAtomOfClut2;
  int nNumAtomOfClutminousone;
  int num,num2;
  double Coord[3],Coord2[3];
  double RotatnNumtonNumMiOn[3][3];
  double mat2[3][3];

  for (alpha=0;alpha<3;++alpha)
    for (alpha2=0;alpha2<3;++alpha2)
      RotatnNumtonNumMiOn[alpha][alpha2] = 0.0;

  for (alpha=0;alpha<3;++alpha)
    for (alpha2=0;alpha2<3;++alpha2)
      for (alpha3=0;alpha3<3;++alpha3)
	RotatnNumtonNumMiOn[alpha][alpha2]
	+= clt[nNumCltminousone].trans_A_to_CN[0][alpha][alpha3]*clt[nNumClt].trans_A_to_CN[0][alpha2][alpha3];

  nNumAtomOfClut = clt[nNumClt].origin_atom_a-1;
  nNumAtomOfClutminousone = clt[nNumCltminousone].origin_atom_a-1;

  for(alpha=0;alpha<3;++alpha)
    Coord[alpha]=( crd[nNumAtomOfClut][alpha] -crd[nNumAtomOfClutminousone][alpha]);

  for(alpha=0;alpha<3;++alpha) Coord2[alpha]=0.0;

  for(alpha=0;alpha<3;++alpha)
    for(alpha2=0;alpha2<3;++alpha2)
      Coord2[alpha]+= clt[nNumCltminousone].trans_A_to_CN[0][alpha][alpha2]*Coord[alpha2];

  for(alpha=0;alpha<3;++alpha)
    for(alpha2=0;alpha2<3;++alpha2)
      clt[nNumClt].TransMatrix[0][alpha][alpha2]=RotatnNumtonNumMiOn[alpha][alpha2];


  for(alpha=3;alpha<6;++alpha)
    for(alpha2=3;alpha2<6;++alpha2)
      clt[nNumClt].TransMatrix[0][alpha][alpha2]
	=RotatnNumtonNumMiOn[alpha-3][alpha2-3];

  for(alpha=3;alpha<6;++alpha)
    for(alpha2=0;alpha2<3;++alpha2)
      clt[nNumClt].TransMatrix[0][alpha][alpha2]=0.0;

  mat2[0][0]= 0.0;
  mat2[0][1]=-Coord2[2];
  mat2[0][2]= Coord2[1];
  mat2[1][0]= Coord2[2];
  mat2[1][1]= 0.0;
  mat2[1][2]=-Coord2[0];
  mat2[2][0]=-Coord2[1];
  mat2[2][1]= Coord2[0];
  mat2[2][2]= 0.0;

  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      clt[nNumClt].TransMatrix[0][i][j+3] = 0.0;

  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	clt[nNumClt].TransMatrix[0][i][j+3] += mat2[i][k]*RotatnNumtonNumMiOn[k][j];

}

