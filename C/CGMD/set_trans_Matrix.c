#include <stdio.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"

//void set_pseduo_trans_matrix(int nNumClut);

//void set_delts_matrix(void);

double Coord2[3];

// 座標系変換行列の作成を行う関数
void set_trans_Matrix(int nNumClt,
					  int nNumClutOrigBranch)
{
	int alpha;
	int alpha2;

	int nNumClutParent;

	nNumClutParent = clust[nNumClt].nNumClutOfParent-1;

//	set_delts_matrix();

	// 0番目の剛体のとき
	if (nNumClt == 0)
	{
		for (alpha=0;alpha</*3*/6;++alpha)
		{
			for (alpha2=0;alpha2</*3*/6;++alpha2)
			{
				clust[nNumClt].TransMatrix[0][alpha][alpha2]= 0.0;
			}
		}
	}
	else
	{
		sub_set_trans_Matrix(nNumClt, nNumClutParent);
	}

	//	set_pseduo_trans_matrix(nNumClt);
}

// 座標系変換行列の作成の補助を行う関数
void sub_set_trans_Matrix(int nNumClt,
	                      int nNumCltminousone)
{
	int alpha,alpha2,alpha3,i,j,k;
	int nNumAtomOfClut,nNumAtomOfClut2;
	int nNumAtomOfClutminousone;
	int num,num2;
	double Coord[3];
	double RotatnNumtonNumMiOn[3][3];
	double mat2[3][3];

	FILE *outtest;

	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			RotatnNumtonNumMiOn[alpha][alpha2] = 0.0;
		}
	}

	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			for (alpha3=0;alpha3<3;++alpha3)
			{
//				// Rot_0_to_N * Rot_n-1_to_0
				RotatnNumtonNumMiOn[alpha][alpha2]
				+=  clust[nNumCltminousone].trans_A_to_CN[0][alpha][alpha3]
			       *clust[nNumClt].trans_A_to_CN[0][alpha2][alpha3];
//				RotatnNumtonNumMiOn[alpha][alpha2]
//				+=  clust[nNumCltminousone].trans_A_to_CN[0][alpha3][alpha]
//			       *clust[nNumClt].trans_A_to_CN[0][alpha3][alpha2];
			}
		}
	}

	nNumAtomOfClut = clust[nNumClt].origin_atom_a-1;
	nNumAtomOfClutminousone = clust[nNumCltminousone].origin_atom_a-1;

	// 座標の取得
	for(alpha=0;alpha<3;++alpha)
	{
		Coord[alpha]/*m*/
		=( prot.coord[nNumAtomOfClut][alpha]
		  -prot.coord[nNumAtomOfClutminousone][alpha])/*A*/
/////////////////////////////////////////////////////////////////////
//		  +prot.coord[nNumAtomOfClutminousone][alpha]/*A*/
/////////////////////////////////////////////////////////////////////
		/**1.0e-10*/;
	}

	// 座標の取得
	for(alpha=0;alpha<3;++alpha)
	{
		Coord2[alpha]=0.0;
	}

	for(alpha=0;alpha<3;++alpha)
	{
		for(alpha2=0;alpha2<3;++alpha2)
		{
//			Coord2[alpha]
//				+=clust[nNumClt].trans_A_to_CN[0][alpha][alpha2]
//			     *Coord[alpha2];
			Coord2[alpha]
				+= clust[nNumCltminousone].trans_A_to_CN[0][alpha][alpha2]
			      *Coord[alpha2];
		}
	}

//	// 原子番号の取得
//	nNumAtomOfClut = clust[nNumClt].origin_atom_a-1;
//	nNumAtomOfClutminousone = clust[nNumCltminousone].origin_atom_a-1;
//	nNumAtomOfClut2 = clust[nNumCltminoustwo].terminal_atom_a[0]-1;
//
//	num = nNumAtomOfClut - nNumAtomOfClut2;
//	num2 = nNumAtomOfClutminousone - nNumAtomOfClut2;
//
//	// 座標の取得
//	for(alpha=0;alpha<3;++alpha)
//	{
//		Coord2[alpha]/*m*/
//		=( clust[nNumCltminousone].xoord_clust/*[0]*/[num2][alpha]
//		  -clust[nNumCltminousone].xoord_clust/*[0]*/[num][alpha])/*A*/
//         /**1.0e-10*/;
//	}

	// 座標系変換行列の左上の作成
	for(alpha=0;alpha<3;++alpha)
	{
		for(alpha2=0;alpha2<3;++alpha2)
		{
			clust[nNumClt].TransMatrix[0][alpha][alpha2]
			         =RotatnNumtonNumMiOn[alpha][alpha2];
//			clust[nNumClt].TransMatrix[0][alpha][alpha2]
//			         =ident[alpha][alpha2];
		}
	}

	// 座標系変換行列の右下の作成
	for(alpha=3;alpha<6;++alpha)
	{
		for(alpha2=3;alpha2<6;++alpha2)
		{
			clust[nNumClt].TransMatrix[0][alpha][alpha2]
			     =RotatnNumtonNumMiOn[alpha-3][alpha2-3];
//			clust[nNumClt].TransMatrix[0][alpha][alpha2]
//			     =ident[alpha-3][alpha2-3];
		}
	}

	// 座標系変換行列の左下の作成
	for(alpha=3;alpha<6;++alpha)
	{
		for(alpha2=0;alpha2<3;++alpha2)
		{
			clust[nNumClt].TransMatrix[0][alpha][alpha2]=0.0;
		}
	}

//	// 座標系変換行列の右上の作成
//	clust[nNumClt].TransMatrix[0][0][3]= 0.0;
//	clust[nNumClt].TransMatrix[0][1][3]=-Coord2[2]/*m*/;
//	clust[nNumClt].TransMatrix[0][2][3]= Coord2[1]/*m*/;
//	clust[nNumClt].TransMatrix[0][0][4]= Coord2[2]/*m*/;
//	clust[nNumClt].TransMatrix[0][1][4]= 0.0;
//	clust[nNumClt].TransMatrix[0][2][4]=-Coord2[0]/*m*/;
//	clust[nNumClt].TransMatrix[0][0][5]=-Coord2[1]/*m*/;
//	clust[nNumClt].TransMatrix[0][1][5]= Coord2[0]/*m*/;
//	clust[nNumClt].TransMatrix[0][2][5]= 0.0;

//	// 座標系変換行列の右上の作成
//	mat2[0][0]= 0.0;
//	mat2[1][0]=-Coord2[2]/*m*/;
//	mat2[2][0]= Coord2[1]/*m*/;
//	mat2[0][1]= Coord2[2]/*m*/;
//	mat2[1][1]= 0.0;
//	mat2[2][1]=-Coord2[0]/*m*/;
//	mat2[0][2]=-Coord2[1]/*m*/;
//	mat2[1][2]= Coord2[0]/*m*/;
//	mat2[2][2]= 0.0;

	// 座標系変換行列の右上の作成
	mat2[0][0]= 0.0;
	mat2[0][1]=-Coord2[2]/*m*/;
	mat2[0][2]= Coord2[1]/*m*/;
	mat2[1][0]= Coord2[2]/*m*/;
	mat2[1][1]= 0.0;
	mat2[1][2]=-Coord2[0]/*m*/;
	mat2[2][0]=-Coord2[1]/*m*/;
	mat2[2][1]= Coord2[0]/*m*/;
	mat2[2][2]= 0.0;

	for (i=0;i<3;++i)
	{
		for (j=0;j<3;++j)
		{
			clust[nNumClt].TransMatrix[0][i][j+3] = 0.0;
		}
	}

	for (i=0;i<3;++i)
	{
		for (j=0;j<3;++j)
		{
			for (k=0;k<3;++k)
			{
				clust[nNumClt].TransMatrix[0][i][j+3] += mat2[i][k]*RotatnNumtonNumMiOn[k][j];
			}
		}
	}

// test
//	clust[2].TransMatrix[0][2][2] = -1.0;
//
}

void set_pseduo_trans_matrix(nNumClut) {
  int alpha,alpha2;
  int nNumAtomOfClut,nNumAtomOfClutminousone;
  int nNumCltminousone;
  double Coord[3],Coord2[3];

  for (alpha=0;alpha<3;++alpha) {
    for (alpha2=0;alpha2<3;++alpha2) {
      clust[nNumClut].PsedoTransMatrix[alpha][alpha2] = clust[nNumClut].TransMatrix[0][alpha][alpha2];
    }
  }

  nNumCltminousone = clust[nNumClut].nNumClutOfParent-1;
  nNumAtomOfClut = clust[nNumClut].origin_atom_a-1;
  nNumAtomOfClutminousone = clust[nNumCltminousone].origin_atom_a-1;

  // 座標の取得
  for(alpha=0;alpha<3;++alpha) {
    Coord[alpha]/*m*/ =(  prot.coord[nNumAtomOfClut][alpha] - prot.coord[nNumAtomOfClutminousone][alpha])/*A*/;
  }
  
  // 座標の取得
  for(alpha=0;alpha<3;++alpha) {
    Coord2[alpha]=0.0;
  }

  for(alpha=0;alpha<3;++alpha) {
    for(alpha2=0;alpha2<3;++alpha2) {
      Coord2[alpha] += clust[nNumCltminousone].trans_A_to_CN[0][alpha][alpha2]*Coord[alpha2];
    }
  }

  clust[nNumClut].PsedoTransMatrix[0][3] = Coord2[0];
  clust[nNumClut].PsedoTransMatrix[1][3] = Coord2[1];
  clust[nNumClut].PsedoTransMatrix[2][3] = Coord2[2];
  
  for (alpha2=0;alpha2<3;++alpha2) {
    clust[nNumClut].PsedoTransMatrix[3][alpha2] = 0.0;
  }

  clust[nNumClut].PsedoTransMatrix[3][alpha2] = 1.0;
}

void set_delts_matrix(void)
{
	int i,j;

	for (i=0;i</*3*/4;++i)
	{
		for (j=0;j</*3*/4;++j)
		{
			delta_matrix[i][j] = 0.0;
		}
	}

	delta_matrix[0][1] = -1.0;
	delta_matrix[1][0] = 1.0;

}
