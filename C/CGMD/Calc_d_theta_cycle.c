#include <stdio.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "force.h"
#include "MD.h"
#include "BD.h"

//BD///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// ���̃X�e�b�v�� 2 �ʊp���v�Z���s���֐�
void calc_d_theta_cycle(void)
{
	int i;
	int nNumClut;

	int nNumClutOrigBranch = 0;

	// ���݂̍\���ł̃^���p�N���̃|�e���V�����G�l���M�[�A
	// �^���p�N���ɋy�ڂ��͂��v�Z
	calc_force(0);

	// ���݂̍\���ł̃^���p�N���� ABI �ABias Force ���v�Z ( ���[ -> �n�_ )
	for (nNumClut=prot.DOF-1; nNumClut>=0; --nNumClut)
	{
		// ���݂̍\���ł̃^���p�N���� ABI ���v�Z
		calc_ABA_cycle(nNumClut);
		// ���݂̍\���ł̃^���p�N���� Bias Force ���v�Z
		calc_Bias_force_cycle(nNumClut);
	}

	// ���� 0 �Ԗڂ� spatial velocity �̐ݒ�
	for(i=0; i<6; ++i)
	{
		clust[0].sp_velo[i] = 0.0;
	}

	// ���݂̍\���ł̃^���p�N���� spatial velocity ���v�Z ( �n�_ -> ���[ )
	for(nNumClut = 0; nNumClut < prot.DOF; ++nNumClut)
	{
		if (clust[nNumClut].num_branch > 1)
		{
			nNumClutOrigBranch = nNumClut;
		}

		// �\���q�����x�̌v�Z���s��
		set_predict_velo(nNumClut, nNumClutOrigBranch);
		// spatial acceleration �̌v�Z���s��
		calc_sp_velo_cycle(nNumClut);
	}
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
