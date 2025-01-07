
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABAb.h"
//#include "EF.h"

void ABAbs_trans_Matrix(CLTb *clt,int nNumClut_all,int num_atom_all,double *crd) {
  int i;
  int nNumClt;
  int alpha,alpha2;  
  int nNumClutParent;
  int *num_index_terminal_atom_a;

  num_index_terminal_atom_a=(int)gcemalloc(sizeof(int)*nNumClut_all);
  for (nNumClt=0;nNumClt<nNumClut_all;++nNumClt) {
    num_index_terminal_atom_a[nNumClt]=0;
  }

  for (nNumClt=0;nNumClt<nNumClut_all;++nNumClt) {
    nNumClutParent = clt[nNumClt].nNumClutOfParent-1;

    if (nNumClt==0)
      for (alpha=0;alpha<6;++alpha)
	for (alpha2=0;alpha2<6;++alpha2)
	  clt[nNumClt].TM[alpha][alpha2]= 0.0;
    else {
      i=num_index_terminal_atom_a[nNumClutParent];
      ++num_index_terminal_atom_a[nNumClutParent];
      sub_set_trans_Matrix(clt[nNumClt].TM,
			   nNumClt,clt[nNumClt].trans_A_to_CN,clt[nNumClt].origin_atom_a,
			   nNumClutParent,clt[nNumClutParent].trans_A_to_CN,
			   clt[nNumClutParent].origin_atom_a/*clt[nNumClutParent].terminal_atom_a[i]*/,
			   crd);
    }
  }
}

void sub_set_trans_Matrix(double TM[6][6],
			  int nNumClt,          double trans_A_to_CN[3][3],          int origin_atom_a,
			  int nNumCltminousone, double trans_A_to_CN_minousone[3][3],int origin_atom_a_minousone,
			  double *crd) {
  int alpha,alpha2,alpha3,i,j,k;
  int nNumAtomOfClut,nNumAtomOfClut2;
  int nNumAtomOfClutminousone;
  int num,num2;
  double Coord[3];
  double Coord2[3];
  double RotatnNumtonNumMiOn[3][3];
  double mat2[3][3];

  for (alpha=0;alpha<3;++alpha)
    for (alpha2=0;alpha2<3;++alpha2)
      RotatnNumtonNumMiOn[alpha][alpha2] = 0.0;

  for (alpha=0;alpha<3;++alpha)
    for (alpha2=0;alpha2<3;++alpha2)
      for (alpha3=0;alpha3<3;++alpha3)
	RotatnNumtonNumMiOn[alpha][alpha2]
	  +=  trans_A_to_CN_minousone[alpha][alpha3]*trans_A_to_CN[alpha2][alpha3];
	
  nNumAtomOfClut = origin_atom_a-1;
  nNumAtomOfClutminousone = origin_atom_a_minousone-1;

  for(alpha=0;alpha<3;++alpha) Coord[alpha]=crd[nNumAtomOfClut*3+alpha]-crd[nNumAtomOfClutminousone*3+alpha];

  for(alpha=0;alpha<3;++alpha) Coord2[alpha]=0.0;

  for(alpha=0;alpha<3;++alpha)
    for(alpha2=0;alpha2<3;++alpha2)
      Coord2[alpha]+=trans_A_to_CN_minousone[alpha][alpha2]*Coord[alpha2];

  for(alpha=0;alpha<3;++alpha)
      for(alpha2=0;alpha2<3;++alpha2)
	TM[alpha][alpha2]=RotatnNumtonNumMiOn[alpha][alpha2];

  for(alpha=3;alpha<6;++alpha)
    for(alpha2=3;alpha2<6;++alpha2)
      TM[alpha][alpha2]=RotatnNumtonNumMiOn[alpha-3][alpha2-3];

  for(alpha=3;alpha<6;++alpha)
    for(alpha2=0;alpha2<3;++alpha2)
      TM[alpha][alpha2]=0.0;

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
      TM[i][j+3] = 0.0;

  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	TM[i][j+3] += mat2[i][k]*RotatnNumtonNumMiOn[k][j];

}
