#include <stdio.h>
#include <stdlib.h>

#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"
#include "MD.h"
#include "physics.h"

// verlet �@�Ŏ��̃X�e�b�v�̊p���x��
// �v�Z���s��
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
