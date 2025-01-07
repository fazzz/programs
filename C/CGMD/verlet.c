#include <stdio.h>
#include <stdlib.h>

#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"
#include "MD.h"
#include "physics.h"

// verlet 法で次のステップの角速度の
// 計算を行う
int verlet(int nNumStep)
{
	int nNumClut;

	calc_d_theta_cycle();

	for (nNumClut=1;nNumClut<prot.DOF;++nNumClut)
	{
		clust[nNumClut].now_deltadihedang[0]
		 = deltat*clust[nNumClut].ddihedang[0];
	}
	return nNumStep;
}
