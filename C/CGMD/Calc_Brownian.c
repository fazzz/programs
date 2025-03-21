#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"
#include "physics.h"
#include "MD.h"
#include "force.h"
#include "BD.h"
#include "genrand.h"

// 揺動の計算を行う関数
void Calc_Brownian(void)
{
	int alpha;
	int nNumClut;


	for (nNumClut = 0; nNumClut < prot.DOF; ++nNumClut)
	{
		// 拡散係数の設定を行う
		set_diffusion_tensor(nNumClut);
		// Brownian 力の計算
	 	Calc_Brownian_force(nNumClut);
	}

}

// 拡散係数の設定を行う関数
void set_diffusion_tensor(int nNumClut)
{
  //	clust[nNumClut].diffusion_tensor_tra = k_B*T_Kelvin_Now
  //					               *clust[nNumClut].friction_tensor_tra;
  //
  //clust[nNumClut].diffusion_tensor_rot = k_B*T_Kelvin_Now
  //					               *clust[nNumClut].friction_tensor_rot;
}

// Brownian 力の計算
void Calc_Brownian_force(int nNumClut)
{
	int i_c;
	int alpha;

	int nNumOrigc;
	int nNumAbsoc;
	double q_c[MAXA][3];

/*	for(i_c=0;i_c < clust[nNumClut].num_atom_clust;++i_c)
	{*/
		// 並進部分の計算
		for (alpha = 0; alpha < 3; ++alpha)
		{
			for (nNumClut = 0; nNumClut < prot.DOF; ++nNumClut)
			{
			  //				clust[nNumClut].f_c.f_Brownian[alpha+3]
			  // = sqrt(clust[nNumClut].diffusion_tensor_tra
			  //       *clust[nNumClut].InertiaMatrix[3][3])
			  //   *nrmd_Box_Muller_method(1, 0, 1);
			}
		}

		nNumOrigc = clust[nNumClut].origin_atom_a;

		nNumAbsoc = clust[nNumClut].origin_atom_a-1 + i_c;

		for(alpha=0;alpha<3;++alpha)
		{
			q_c[i_c][alpha]
			 = q[nNumOrigc][nNumAbsoc][alpha]/*A*/*1.0e-10/*m/A*/;/*m*/
		}


///////////////////////////////////////////////////////////////////////////////
		// 回転部分の計算
		for (alpha = 0; alpha < 3; ++alpha)
		{
			for (nNumClut = 0; nNumClut < prot.DOF; ++nNumClut)
			{
				clust[nNumClut].f_c.f_Brownian[alpha]
				 = clust[nNumClut].f_c.f_Brownian[alpha+3];
			}
		}
/*	}*/
///////////////////////////////////////////////////////////////////////////////
}
