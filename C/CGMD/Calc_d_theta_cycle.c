#include <stdio.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "force.h"
#include "MD.h"
#include "BD.h"

//BD///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// 次のステップの 2 面角を計算を行う関数
void calc_d_theta_cycle(void)
{
	int i;
	int nNumClut;

	int nNumClutOrigBranch = 0;

	// 現在の構造でのタンパク質のポテンシャルエネルギー、
	// タンパク質に及ぼす力を計算
	calc_force(0);

	// 現在の構造でのタンパク質の ABI 、Bias Force を計算 ( 末端 -> 始点 )
	for (nNumClut=prot.DOF-1; nNumClut>=0; --nNumClut)
	{
		// 現在の構造でのタンパク質の ABI を計算
		calc_ABA_cycle(nNumClut);
		// 現在の構造でのタンパク質の Bias Force を計算
		calc_Bias_force_cycle(nNumClut);
	}

	// 剛体 0 番目の spatial velocity の設定
	for(i=0; i<6; ++i)
	{
		clust[0].sp_velo[i] = 0.0;
	}

	// 現在の構造でのタンパク質の spatial velocity を計算 ( 始点 -> 末端 )
	for(nNumClut = 0; nNumClut < prot.DOF; ++nNumClut)
	{
		if (clust[nNumClut].num_branch > 1)
		{
			nNumClutOrigBranch = nNumClut;
		}

		// 予測子加速度の計算を行う
		set_predict_velo(nNumClut, nNumClutOrigBranch);
		// spatial acceleration の計算を行う
		calc_sp_velo_cycle(nNumClut);
	}
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
