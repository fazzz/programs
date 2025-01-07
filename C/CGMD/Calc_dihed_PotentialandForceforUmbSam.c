#include <stdio.h>
#include <math.h>
#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"
#include "MD.h"
#include "force.h"
#include "UmsSan.h"



// 2 面角相互作用の計算を行う関数
void Calc_dihed_Umbrella_Potential(void)
{
	int nNumClut,nNumDihed;

	double dihedang;
	double P;
	double N;

	nNumClut = atom_pair_US[5]-1;
	dihedang = pick_dihed_one_clust(atom_pair_US[0]-1,
					atom_pair_US[1]-1,
					atom_pair_US[2]-1,
					atom_pair_US[3]-1,
					nNumDihed);

	// ポテンシャルエネルギーの計算
	P = V_K_US*(dihedang-dihed_ref_US)*(dihedang-dihed_ref_US);

	// 力の計算
	N = -2.0*V_K_US*(dihedang-dihed_ref_US)/1.660539e-27/1.0e-20*1.0e-24/2.3889e-4;

	potential_US = P;
	clust[nNumClut].f_c.f_dihed += N;
}

