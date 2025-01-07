#include <stdio.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "force.h"
#include "MD.h"
#include "BD.h"

// BD /////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// 現在の構造でのタンパク質の
// spatial velocity の計算を行う関数
void calc_sp_velo_cycle(int nNumClut)
{
	int i;

	// spatial velocity の初期化を行う
	for(i=0; i<6; ++i)
	{
		clust[nNumClut].sp_velo[i] = 0.0;
	}

	// spatial velocity の計算を行う
	sub_calc_sp_velo_cycle(nNumClut);
}

// 現在の構造でのタンパク質の
// spatial velocity の計算の補助を行う関数
void sub_calc_sp_velo_cycle(int nNumClut)
{
	int i;

	int alpha,alpha2;

	int nNumFreedom;

	double q[3];

	double OmegaOnThisCoord[3];
	double OmegaOnLaboCoord[3];

	// 剛体のspatial velocity の初期化を行う
	for(i=0; i<6; ++i)
	{
		clust[nNumClut].sp_velo[i] = 0.0;
	}


	// 剛体の2面角加速度の初期化を行う
	clust[nNumClut].ddihedang[0] = 0.0;

////////////////////////////////////////////////////////////////////////
	// 剛体の2面角速度の計算を行う_1
	for(i=0; i<6; ++i)
	{
//		clust[nNumClut].ddihedang[0]/*rad/s^2*/
//						-=   ABI[nNumClut].Kalman_Gain_Transpose[i]
//					    	*clust[nNumClut].predict_velo[i]/*rad/s^2*/;
	}

	// 剛体の2面角速度の計算を行う_2
	//0410	clust[nNumClut].ddihedang[0]/*rad/s^2*/ 
	//				+=  zzz[nNumClut].nyu
	//				   /clust[nNumClut].friction_tensor_tra/*rad/s^2*/;
////////////////////////////////////////////////////////////////////////

	// spatial velocity の計算を行う_1
	for(i=0; i<6; ++i)
	{
		clust[nNumClut].sp_velo[i]/*rad/s^2*/
		                 +=   clust[nNumClut].predict_velo[i]/*rad/s^2*/;
	}

/*	// クラスタ座標系での角加速度と
	// 実験室座標系での角加速度の初期化
	for (alpha=0;alpha<3;++alpha)
	{
		OmegaOfThisCoord[alpha] = 0.0;
		OmegaOfAbsoCoord[alpha] = 0.0;
	}

	// このクラスタの自由度の取得
	nNumFreedom = ABI[nNumClut].hingmat-1;

	// クラスタ座標系での角加速度の取得
	OmegaOfThisCoord[nNumFreedom] = clust[nNumClut].ddihedang[0];

	// 実験室座標系での角加速度の計算
	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			OmegaOfAbsoCoord[alpha]
				+= clust[nNumClut].trans_A_to_CN[0][alpha2][alpha]
				  *OmegaOfThisCoord[alpha2];
		}
	}
*/

	// spatial acceleration の計算を行う_2
/*	for (alpha=0;alpha<3;++alpha)
	{
		clust[nNumClut].sp_velo[alpha]/*rad/s^2*/ 
/*		                 += OmegaOfAbsoCoord[alpha];/*rad/s^2*/
/*	}*/
 	// 自由度の選択
	nNumFreedom = ABI[nNumClut].hingmat-1;

 	// クラスタ座標系での角速度の初期化
	for (alpha=0;alpha<3;++alpha)
	{
		OmegaOnThisCoord[alpha] = 0.0;
		OmegaOnLaboCoord[alpha] = 0.0;
	}

 	// クラスタ座標系での角速度の取得
	OmegaOnThisCoord[nNumFreedom] 
	                 = clust[nNumClut].ddihedang[0];

 	// 実験室座標系での角速度
/*	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			OmegaOnLaboCoord[alpha]
					 += clust[nNumClut].trans_A_to_CN[0][alpha2][alpha]
		               *OmegaOnThisCoord[alpha2];
		}
	}*/

	 // 実験室座標系での角速度
	OmegaOnLaboCoord[0]
		=   clust[nNumClut].trans_A_to_CN[0][1][1]
		   *clust[nNumClut].trans_A_to_CN[0][2][0]
		   *clust[nNumClut].ddihedang[0]
		  - clust[nNumClut].trans_A_to_CN[0][1][0]
		   *clust[nNumClut].trans_A_to_CN[0][2][1]
		   *clust[nNumClut].ddihedang[0];

	OmegaOnLaboCoord[1]
		=   clust[nNumClut].trans_A_to_CN[0][0][1]
		   *clust[nNumClut].trans_A_to_CN[0][2][0]
		   *clust[nNumClut].ddihedang[0]
		  - clust[nNumClut].trans_A_to_CN[0][0][0]
		   *clust[nNumClut].trans_A_to_CN[0][2][1]
		   *clust[nNumClut].ddihedang[0];

	OmegaOnLaboCoord[2]
		=   clust[nNumClut].trans_A_to_CN[0][0][0]
		   *clust[nNumClut].trans_A_to_CN[0][1][1]
		   *clust[nNumClut].ddihedang[0]
		  - clust[nNumClut].trans_A_to_CN[0][0][1]
		   *clust[nNumClut].trans_A_to_CN[0][1][0]
		   *clust[nNumClut].ddihedang[0];

	// spatial velocity の設定を行う_2
	for (alpha=0;alpha<3;++alpha)
	{
		clust[nNumClut].sp_velo[alpha] += OmegaOnLaboCoord[alpha];
	}
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
