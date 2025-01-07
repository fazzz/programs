#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "physics.h"

// 予測加速度の計算の補助を行う関数_2
void sub_sub_set_predict_acc(int nNumClt, int nNumCltminusone);

// 予測加速度の計算を行う関数
void set_predict_acc(int nNumClt, int nNumClutOrigBranch)
{
	int i;
	int nNumClutOfParent;

	nNumClutOfParent = clust[nNumClt].nNumClutOfParent-1;

	// クラスタが末端のとき
	if (nNumClutOfParent == 0)
	{
		for(i=0;i<6;++i)
		{
			clust[nNumClt].predict_alpha[i] = 0.0;
		}
	}
	else
	{
		sub_set_predict_acc(nNumClt, nNumClutOfParent);
	}
}

// 予測加速度の計算の補助を行う関数_1
void sub_set_predict_acc(int nNumClt, int nNumCltminusone)
{
	int i,j;

	// 剛体の予測加速度の初期化を行う
	for(i=0; i<6; ++i)
	{
		clust[nNumClt].predict_alpha[i] = 0.0;
	}

	// 変換行列の転置行列を乗じる
	for (i=0;i<6;++i)
	{
		for (j=0;j<6;++j)
		{
			clust[nNumClt].predict_alpha[i]
			 +=   clust[nNumClt].TransMatrix[0][j][i]
				* clust[nNumCltminusone].sp_acc[j];
		}
	}
}

// 予測加速度の計算の補助を行う関数_2
//void sub_sub_set_predict_acc(int nNumClt, int nNumCltminusone)
//{
//	int i,j;
//
//	double TransMatrix[6][6];
//
//	// 剛体の予測加速度の初期化を行う
//	for(i=0; i<6; ++i)
//	{
//		clust[nNumClt].predict_alpha[i] = 0.0;
//	}
//
//	for (i=0;i<6;++i)
//	{
//		for (j=0;j<6;++j)
//		{
//			TransMatrix[i][j] = 0.0;
//		}
//	}
//
//	for (i=0;i<6;++i)
//	{
//		for (j=0;j<6;++j)
//		{
//			TransMatrix[i][j] = -1*clust[nNumClt].TransMatrix[0][i][j];
//		}
//	}
//
//	for (i=0;i<6;++i)
//	{
//		TransMatrix[i][i] += 2;
//	}
//
//	// 変換行列の転置行列を乗じる
//	for (i=0;i<6;++i)
//	{
//		for (j=0;j<6;++j)
//		{
//			clust[nNumClt].predict_alpha[i]
//			 +=   TransMatrix[j][i]
//				* clust[nNumCltminusone].sp_acc[j];
//		}
//	}
//}
