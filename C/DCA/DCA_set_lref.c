
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "DCA.h"
#include "EF.h"

void DCAs_local_reference(CLT *clt,int nNumClut_all,int num_atom_all,double *crd) {
  int i,j;
  int nNumClut;
  int nNumClutOrigBranch;
  int nNumClutOfParent;
  int nNumClut2;
  int num,nNumClutdummy;
  int flag;

  for (nNumClut=0;nNumClut<nNumClut_all;++nNumClut) {
    nNumClutOfParent = clt[nNumClut].nNumClutOfParent-1;

    if(clt[nNumClut].num_branch > 1) nNumClutOrigBranch = nNumClut;

    if (nNumClut == 0)  {
      clt[nNumClut].xoord=(double *)gcemalloc(sizeof(double)*num_atom_all*3);
      sub_trans_A_to_CN(clt[nNumClut].trans_A_to_CN,clt[nNumClut].xoord,
			0,clt[nNumClut].origin_atom_a,
			-1,0,
			clt[nNumClutOfParent].origin_atom_a,
			num_atom_all,num_atom_all,crd);
    }
    else if(clt[nNumClut].join > 0) {
      num=0;
      for (nNumClutdummy=nNumClut;clt[nNumClutdummy].join!=clt[nNumClut].join-1;++nNumClutdummy)
	num+=clt[nNumClutdummy].num_atom_clust;
      clt[nNumClut].xoord=(double *)gcemalloc(sizeof(double)*num_atom_all*3);
      sub_trans_A_to_CN(clt[nNumClut].trans_A_to_CN,
			clt[nNumClut].xoord,
			nNumClut,clt[nNumClut].origin_atom_a,
			clt[nNumClutOfParent].origin_atom_a,
			nNumClutOfParent,clt[nNumClutOfParent].origin_atom_a,
			num,num_atom_all,crd);
    }
    else {
      clt[nNumClut].xoord=(double *)gcemalloc(sizeof(double)*num_atom_all*3);
      sub_trans_A_to_CN(clt[nNumClut].trans_A_to_CN,
			clt[nNumClut].xoord,
			nNumClut,clt[nNumClut].origin_atom_a,
			nNumClutOfParent,clt[nNumClutOfParent].origin_atom_a,
			clt[nNumClutOfParent].origin_atom_a,
			num_atom_all-clt[nNumClut].origin_atom_a+1,
			num_atom_all,crd);
    }
  }
}
  
void sub_trans_A_to_CN(double Mtrans_A_to_CN[3][3], double *crd_local,
		       int nNumCltTar,   int origin_atom_a,
		       int nNumCltCoo,   int terminal_atom_a,
		       int xy_set_atom_a,
		       int nNumAtom,     int num_atom,
		       double *crd){
  int i;  
  int alpha; 
  int nNmAtomOfCN_A;
  int nNmAtomOfCN_1_A;
  
  double CN_A[3];
  double HN_A[3];
  double CN_1_A[3];
  
  double *mat;

  double ii[3];
  double jj[3];
  double kk[3];
  
  double jj_x;
  double jj_y;
  double jj_z;
  
  double DisOfCN_1_CN=0.0;
  double DisOfHN_CN=0.0;
  double SnCN_1_CN_HN=0.0;
  double CsCN_1_CN_HN=0.0;
  
  mat=(double *)gcemalloc(sizeof(double)*num_atom*3);

  nNmAtomOfCN_A = origin_atom_a-1;
  if (nNumCltCoo != -1)  nNmAtomOfCN_1_A = terminal_atom_a-1;
  else nNmAtomOfCN_1_A = 0;
  
  for(alpha=0;alpha<3;++alpha)	{
    if (nNumCltCoo != -1)	{
      CN_A[alpha]=crd[nNmAtomOfCN_A*3+alpha];
      HN_A[alpha]=crd[(nNmAtomOfCN_A+1)*3+alpha];
      //      HN_A[alpha]=crd[(xy_set_atom_a-1)*3+alpha];
    }
    else {
      CN_A[alpha]=crd[nNmAtomOfCN_A*3+alpha];
      HN_A[alpha]=crd[2*3+alpha];
    }
    if (nNumCltCoo != -1)	{
      CN_1_A[alpha]=crd[nNmAtomOfCN_1_A*3+alpha];
    }
    else{
      CN_1_A[alpha]=0.0;
    }
  }
  
  for(alpha=0;alpha<3;++alpha)	{
    kk[alpha] = CN_1_A[alpha]-CN_A[alpha];
    DisOfCN_1_CN += (CN_1_A[alpha]-CN_A[alpha])*(CN_1_A[alpha]-CN_A[alpha]);
    DisOfHN_CN += (HN_A[alpha]-CN_A[alpha])*(HN_A[alpha]-CN_A[alpha]);
    CsCN_1_CN_HN += (CN_1_A[alpha]-CN_A[alpha])*(HN_A[alpha]-CN_A[alpha]);
  }
  
  DisOfCN_1_CN = sqrt(DisOfCN_1_CN);
  DisOfHN_CN = sqrt(DisOfHN_CN);
  CsCN_1_CN_HN = CsCN_1_CN_HN/(DisOfCN_1_CN*DisOfHN_CN);
  SnCN_1_CN_HN = 1.0-CsCN_1_CN_HN*CsCN_1_CN_HN;
  SnCN_1_CN_HN = sqrt(SnCN_1_CN_HN);
  
  for(alpha=0;alpha<3;++alpha) kk[alpha]=kk[alpha]/DisOfCN_1_CN;
  
  jj[0]=((CN_1_A[1]-CN_A[1])*(HN_A[2]-CN_A[2])
	 -(CN_1_A[2]-CN_A[2])*(HN_A[1]-CN_A[1]))
    /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);
  jj[1]=((CN_1_A[2]-CN_A[2])*(HN_A[0]-CN_A[0])
	 -(CN_1_A[0]-CN_A[0])*(HN_A[2]-CN_A[2]))
    /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);
  jj[2]=((CN_1_A[0]-CN_A[0])*(HN_A[1]-CN_A[1])
	 -(CN_1_A[1]-CN_A[1])*(HN_A[0]-CN_A[0]))
    /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);
  
  ii[0]=jj[1]*kk[2]-jj[2]*kk[1];
  ii[1]=jj[2]*kk[0]-jj[0]*kk[2];
  ii[2]=jj[0]*kk[1]-jj[1]*kk[0];
  
  for(i=0; i<nNumAtom; ++i){
    for(alpha=0;alpha<3;++alpha){
      if (nNumCltCoo == -1)
	mat[i*3+alpha] = crd[i*3+alpha]-CN_A[alpha];
      else
	mat[i*3+alpha] = crd[(nNmAtomOfCN_A+i)*3+alpha]-CN_A[alpha];	    
    }
  }
  
  for(alpha=0;alpha<3;++alpha){
    Mtrans_A_to_CN[0][alpha]=ii[alpha];
    Mtrans_A_to_CN[1][alpha]=jj[alpha];
    Mtrans_A_to_CN[2][alpha]=kk[alpha];
  }
  
  for(i=0; i < nNumAtom; ++i){
    crd_local[i*3+0]
      =   Mtrans_A_to_CN[0][0]*mat[i*3+0]
      + Mtrans_A_to_CN[0][1]*mat[i*3+1]
      + Mtrans_A_to_CN[0][2]*mat[i*3+2];
    
    crd_local[i*3+1]
      =   Mtrans_A_to_CN[1][0]*mat[i*3+0]
      + Mtrans_A_to_CN[1][1]*mat[i*3+1]
      + Mtrans_A_to_CN[1][2]*mat[i*3+2];
    
    crd_local[i*3+2]
      =   Mtrans_A_to_CN[2][0]*mat[i*3+0]
      + Mtrans_A_to_CN[2][1]*mat[i*3+1]
      + Mtrans_A_to_CN[2][2]*mat[i*3+2];
  }
}
