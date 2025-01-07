#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "physics.h"

// 予測速度の計算の補助を行う関数_1
void sub_set_predict_velo(int nNumClt, int nNumCltminusone);
// 予測速度の計算の補助を行う関数_2
void sub_sub_set_predict_velo(int nNumClt, int nNumCltminusone);

//BD///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// 予測速度の計算を行う関数
void set_predict_velo(int nNumClt, int nNumClutOrigBranch)
{
	int i;

	// 0番目の剛体のとき
	if (nNumClt == 0)
	{
		for (i=0;i<6;++i)
		{
			clust[nNumClt].predict_velo[i] = 0.0;
		}
	}
	// 終端の剛体の予測速度の計算を行う
	else if(clust[nNumClt].terminal == TERMINAL)
	{
		sub_sub_set_predict_velo(nNumClt, nNumClutOrigBranch);
	}
	// 分岐の剛体の予測速度の計算を行う
	else if (clust[nNumClt-1].terminal == TERMINAL)
	{
		sub_set_predict_velo(nNumClt, nNumClutOrigBranch);
	}
	// 通常の剛体の予測速度の計算を行う
	else
	{
		sub_set_predict_velo(nNumClt, nNumClt-1);
	}
}

// 予測速度の計算の補助を行う関数_1
void sub_set_predict_velo(int nNumClt, int nNumCltminusone)
{
	int i,j;

	// 剛体の予測速度の初期化を行う
	for(i=0; i<6; ++i)
	{
		clust[nNumClt].predict_velo[i] = 0.0;
	}

	// 変換行列の転置行列を乗じる
	for (i=0;i<6;++i)
	{
		for (j=0;j<6;++j)
		{
			clust[nNumClt].predict_velo[i]
			 +=   clust[nNumClt].TransMatrix[0][j][i]
				* clust[nNumCltminusone].sp_velo[j];
		}
	}
}

// 予測加速度の計算の補助を行う関数_2
void sub_sub_set_predict_velo(int nNumClt, int nNumCltminusone)
{
	int i,j;

	double TransMatrix[6][6];

	// 剛体の予測速度の初期化を行う
	for(i=0; i<6; ++i)
	{
		clust[nNumClt].predict_velo[i] = 0.0;
	}

	for (i=0;i<6;++i)
	{
		for (j=0;j<6;++j)
		{
			TransMatrix[i][j] = 0.0;
		}
	}

	for (i=0;i<6;++i)
	{
		for (j=0;j<6;++j)
		{
			TransMatrix[i][j] = -1.0*clust[nNumClt].TransMatrix[0][i][j];
		}
	}

	for (i=0;i<6;++i)
	{
		TransMatrix[i][i] += 2.0;
	}

	// 変換行列の転置行列を乗じる
	for (i=0;i<6;++i)
	{
		for (j=0;j<6;++j)
		{
			clust[nNumClt].predict_velo[i]
			 +=   TransMatrix[j][i]
				* clust[nNumCltminusone].sp_velo[j];
		}
	}
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
