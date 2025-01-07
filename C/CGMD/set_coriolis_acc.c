#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "physics.h"
#include "MD.h"

// corioli 加速度の計算を行う関数
void set_coriolis_acc(int nNumClut  // 剛体のインデックス
                     ) {
  int i;
  int nNumClutOfParent;

  nNumClutOfParent = clust[nNumClut].nNumClutOfParent-1;

//	// クラスタが末端のとき
//	if (nNumClutOfParent == 0)
//	{
//		for(i=0;i<6;++i)
//		{
//			clust[nNumClut].Coriolis_acc[i] = 0.0;
//		}
//	}
  // クラスタが末端のとき
  if (nNumClut == 0){
    if (TermMoveMode==OFF) {
      for(i=0;i<6;++i) {
	clust[nNumClut].Coriolis_acc[i] = 0.0;
      }
    }
    else {
      sub_set_coriolis_acc(nNumClut, nNumClutOfParent, -1);
      for(i=0;i<6;++i) {
	clust[nNumClut].Coriolis_acc[i] = 0.0;
      }
    }
  }
  else {
    sub_set_coriolis_acc(nNumClut, nNumClutOfParent, 0);
  }
  
}

// corioli 加速度の計算の補助を行う関数
void sub_set_coriolis_acc(int nNumClut,
			  int nNumClutminusone,
			  int nNumBod)
{
	int i,j;
	int alpha;
	int alpha2;
	int alpha3;
	int nNumAtomOfClut;
	int nNumAtomOfClutminousone;

	// 剛体系nNumClutminusone->剛体系nNumClutへの回転マト
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
	double sp_velo_minous_one[3];
	double Coriolis_acc_dummy_dummy[3];

	// corioli 加速度の初期化
	for (i=0;i<6;++i)
	{
		clust[nNumClut].Coriolis_acc[i] = 0.0;
	}

	for (alpha=0;alpha<3;++alpha)
	{
		Omega[alpha] = 0.0;
		vector[alpha] = 0.0;
		vector2[alpha] = 0.0;
		Coriolis_acc_dummy[alpha] = 0.0;
		Coriolis_acc_dummy_up[alpha] = 0.0;
		Coriolis_acc_dummy_up2[alpha] = 0.0;
	}

	for (alpha=0;alpha<3;++alpha)
	{
//		for (alpha2=0;alpha2<3;++alpha2)
//		{
			Omega[alpha] /*+*/= //clust[nNumClut].TransMatrix[0]/*[alpha][alpha2]*/[alpha2][alpha]
				           /***/clust[nNumClut/*minusone*/].sp_velo[alpha/*2*/];
//		}
	}

//	// 上半分(回転部分)の設定
//	Coriolis_acc_dummy_up[0] = /*-*/Omega[1]
//	                                  *clust[nNumClut].ddihedang[0];
//
//	Coriolis_acc_dummy_up[1] = -/*-*/Omega[0]
//	                                  *clust[nNumClut].ddihedang[0];
//
//	Coriolis_acc_dummy_up[2] = 0.0;
//
//	for (i=0;i<3;++i)
//	{
//		for (j=0;j<3;++j)
//		{
//			clust[nNumClut].Coriolis_acc[i] 
//			+= clust[nNumClut].TransMatrix[0][j][i]*Coriolis_acc_dummy_up[j];
//		}
//	}

	ez[0] = 0.0;
	ez[1] = 0.0;
	ez[2] = 1.0;

	for (i=0;i<3;++i)
	{
		for (j=0;j<3;++j)
		{
			Coriolis_acc_dummy_up[i] 
			/*+*/= /*clust[nNumClut].TransMatrix[0][j][i]**/Omega[/*j*/i];
		}
	}
	omegaproduct2[0][0]= 0.0;
	omegaproduct2[0][1]=-Coriolis_acc_dummy_up[2]/*m*/;
	omegaproduct2[0][2]= Coriolis_acc_dummy_up[1]/*m*/;
	omegaproduct2[1][0]= Coriolis_acc_dummy_up[2]/*m*/;
	omegaproduct2[1][1]= 0.0;
	omegaproduct2[1][2]=-Coriolis_acc_dummy_up[0]/*m*/;
	omegaproduct2[2][0]=-Coriolis_acc_dummy_up[1]/*m*/;
	omegaproduct2[2][1]= Coriolis_acc_dummy_up[0]/*m*/;
	omegaproduct2[2][2]= 0.0;

	for (i=0;i<3;++i)
	{
		for (j=0;j<3;++j)
		{
			Coriolis_acc_dummy_up2[i]
				+= omegaproduct2[i][j]*ez[j];
		}
	}

	for (i=0;i<3;++i)
	{
		clust[nNumClut].Coriolis_acc[i] = Coriolis_acc_dummy_up2[i]*clust[nNumClut].ddihedang[0];
	}

//	// 下半分(並進部分)の設定
//	clust[nNumClut].Coriolis_acc[3] =( -clust[nNumClutminusone].sp_velo[1]
//		                               *clust[nNumClutminusone].sp_velo[1]
//									   -clust[nNumClutminusone].sp_velo[2]
//		                               *clust[nNumClutminusone].sp_velo[2] )
//									  *clust[nNumClut].TransMatrix[0][1][5]
//									  + clust[nNumClutminusone].sp_velo[0]
//	                                   *clust[nNumClutminusone].sp_velo[1]
//									  *clust[nNumClut].TransMatrix[0][2][3]
//									  + clust[nNumClutminusone].sp_velo[0]
//	                                   *clust[nNumClutminusone].sp_velo[2]
//									  *clust[nNumClut].TransMatrix[0][0][4];
//
//	clust[nNumClut].Coriolis_acc[4] =  clust[nNumClutminusone].sp_velo[0]
//	                                  *clust[nNumClutminusone].sp_velo[1]
//									 *clust[nNumClut].TransMatrix[0][1][5]
//									 + ( -clust[nNumClutminusone].sp_velo[2]
//									 	 *clust[nNumClutminusone].sp_velo[2]
//									     -clust[nNumClutminusone].sp_velo[0]
//									 	*clust[nNumClutminusone].sp_velo[0] )
//									 *clust[nNumClut].TransMatrix[0][2][3]
//									 + clust[nNumClutminusone].sp_velo[1]
//	                                  *clust[nNumClutminusone].sp_velo[2]
//									  *clust[nNumClut].TransMatrix[0][0][4];
//
//	clust[nNumClut].Coriolis_acc[5] = clust[nNumClutminusone].sp_velo[0]
//	                                 *clust[nNumClutminusone].sp_velo[2]
//									 *clust[nNumClut].TransMatrix[0][1][5]
//									+ clust[nNumClutminusone].sp_velo[2]
//	                                 *clust[nNumClutminusone].sp_velo[1]
//									 *clust[nNumClut].TransMatrix[0][2][3]
//									+ ( -clust[nNumClutminusone].sp_velo[1]
//									    *clust[nNumClutminusone].sp_velo[1]
//									    -clust[nNumClutminusone].sp_velo[0]
//										*clust[nNumClutminusone].sp_velo[0] )
//									  *clust[nNumClut].TransMatrix[0][0][4];

	// 原子番号の取得
	nNumAtomOfClut = clust[nNumClut].origin_atom_a-1;
	if (nNumClutminusone==-1) {
	  vector[alpha] = prot.coord[nNumAtomOfClut][alpha];
	}
	else {
	  nNumAtomOfClutminousone = clust[nNumClutminusone].origin_atom_a-1;
	  for (alpha=0;alpha<3;++alpha){
	    vector[alpha] = prot.coord[nNumAtomOfClut][alpha]
	      -prot.coord[nNumAtomOfClutminousone][alpha];
	  }
	}

	for (alpha=0;alpha<3;++alpha)
	{
		vector2[alpha] = 0.0;
	}

	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			vector2[alpha] += clust[nNumClutminusone].trans_A_to_CN[0]/*[alpha2][alpha]*/[alpha][alpha2]
							 *vector/*[alpha]*/[alpha2]/**1.0e-10*/;
		}
	}

//	// 下半分(並進部分)の設定
//	Coriolis_acc_dummy[0] =( 
//					   -clust[nNumClutminusone].sp_velo[1]
//	                   *clust[nNumClutminusone].sp_velo[1]
//					   -clust[nNumClutminusone].sp_velo[2]
//	                   *clust[nNumClutminusone].sp_velo[2] )
//					  *vector2[0]
//					  + clust[nNumClutminusone].sp_velo[0]
//	                   *clust[nNumClutminusone].sp_velo[1]
//					  *vector2[1]
//					  + clust[nNumClutminusone].sp_velo[0]
//	                   *clust[nNumClutminusone].sp_velo[2]
//					  *vector2[2];
//
//	Coriolis_acc_dummy[1] =  clust[nNumClutminusone].sp_velo[0]
//	                  *clust[nNumClutminusone].sp_velo[1]
//					  *vector2[0]
//					 + ( -clust[nNumClutminusone].sp_velo[2]
//					 	 *clust[nNumClutminusone].sp_velo[2]
//					     -clust[nNumClutminusone].sp_velo[0]
//					 	*clust[nNumClutminusone].sp_velo[0] )
//					  *vector2[1]
//					 + clust[nNumClutminusone].sp_velo[1]
//	                  *clust[nNumClutminusone].sp_velo[2]
//					  *vector2[2];
//
//	Coriolis_acc_dummy[2] = clust[nNumClutminusone].sp_velo[0]
//	                 *clust[nNumClutminusone].sp_velo[2]
//					 *vector2[0]
//					+ clust[nNumClutminusone].sp_velo[2]
//	                 *clust[nNumClutminusone].sp_velo[1]
//					 *vector2[1]
//					+ ( -clust[nNumClutminusone].sp_velo[1]
//					    *clust[nNumClutminusone].sp_velo[1]
//					    -clust[nNumClutminusone].sp_velo[0]
//						*clust[nNumClutminusone].sp_velo[0] )
//					  *vector2[2];


//	omegaproduct[0][0]= 0.0;
//	omegaproduct[0][1]=-clust[nNumClut].sp_velo[2]/*m*/;
//	omegaproduct[0][2]= clust[nNumClut].sp_velo[1]/*m*/;
//	omegaproduct[1][0]= clust[nNumClut].sp_velo[2]/*m*/;
//	omegaproduct[1][1]= 0.0;
//	omegaproduct[1][2]=-clust[nNumClut].sp_velo[0]/*m*/;
//	omegaproduct[2][0]=-clust[nNumClut].sp_velo[1]/*m*/;
//	omegaproduct[2][1]= clust[nNumClut].sp_velo[0]/*m*/;
//	omegaproduct[2][2]= 0.0;

	omegaproduct[0][0]= 0.0;
	omegaproduct[0][1]=-clust[nNumClutminusone].sp_velo[2]/*m*/;
	omegaproduct[0][2]= clust[nNumClutminusone].sp_velo[1]/*m*/;
	omegaproduct[1][0]= clust[nNumClutminusone].sp_velo[2]/*m*/;
	omegaproduct[1][1]= 0.0;
	omegaproduct[1][2]=-clust[nNumClutminusone].sp_velo[0]/*m*/;
	omegaproduct[2][0]=-clust[nNumClutminusone].sp_velo[1]/*m*/;
	omegaproduct[2][1]= clust[nNumClutminusone].sp_velo[0]/*m*/;
	omegaproduct[2][2]= 0.0;
//	omegaproduct[0][0]= 0.0;
//	omegaproduct[0][1]=-clust[nNumClutminusone].ddihedang[0]/*m*/;
//	omegaproduct[0][2]= 0.0/*m*/;
//	omegaproduct[1][0]= clust[nNumClutminusone].ddihedang[0]/*m*/;
//	omegaproduct[1][1]= 0.0;
//	omegaproduct[1][2]= 0.0;/*m*/;
//	omegaproduct[2][0]= 0.0;/*m*/;
//	omegaproduct[2][1]= 0.0;/*m*/;
//	omegaproduct[2][2]= 0.0;


	for (alpha=0;alpha<3;++alpha)
	{
		Coriolis_acc_dummy[alpha] = 0.0;
		Coriolis_acc_dummy_dummy[alpha] = 0.0;
	}

	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			Coriolis_acc_dummy_dummy[alpha] += omegaproduct[alpha][alpha2]*vector2[alpha2];
		}
	}

	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			Coriolis_acc_dummy[alpha] += omegaproduct[alpha][alpha2]*Coriolis_acc_dummy_dummy[alpha2];
		}
	}

//	for (alpha=0;alpha<3;++alpha)
//	{
//		for (alpha2=0;alpha2<3;++alpha2)
//		{
//			matOmegaOmega[alpha][alpha2] = 0.0;
//		}
//	}
//
//	for (alpha=0;alpha<3;++alpha)
//	{
//		for (alpha2=0;alpha2<3;++alpha2)
//		{
//			for (alpha3=0;alpha3<3;++alpha3)
//			{
//				matOmegaOmega[alpha][alpha2] += omegaproduct[alpha][alpha3]
//											   *omegaproduct[alpha3][alpha2];
//			}
//		}
//	}
//
//	for (alpha=0;alpha<3;++alpha)
//	{
//		Coriolis_acc_dummy[alpha] = 0.0;
//	}
//
//	for (alpha=0;alpha<3;++alpha)
//	{
//		for (alpha2=0;alpha2<3;++alpha2)
//		{
//			Coriolis_acc_dummy[alpha] += matOmegaOmega[alpha][alpha2]
//										*vector2[alpha2];
//		}
//	}
//
	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			clust[nNumClut].Coriolis_acc[alpha+3]
							  += clust[nNumClut].TransMatrix[0]/*[alpha][alpha2]*/[alpha2][alpha]
				              *Coriolis_acc_dummy[alpha2];
		}
	}

/////////////////////////////////////////////////////////////////////
	for (i=0;i<3;++i)
	{
		sp_velo_minous_one[i] = 0.0;
	}
	for (i=0;i<3;++i)
	{
		for (j=0;j<3;++j)
		{
			sp_velo_minous_one[i] 
			+= clust[nNumClut].TransMatrix[0][j][i]*clust[nNumClutminusone].sp_velo[j+3];
		}
	}

	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
//			clust[nNumClut].Coriolis_acc[alpha+3]
//								+= omegaproduct[alpha][alpha2]
//								  *(clust[nNumClut].sp_velo[alpha2+3]-sp_velo_minous_one[alpha2]);
		}
	}
/////////////////////////////////////////////////////////////////////
}
