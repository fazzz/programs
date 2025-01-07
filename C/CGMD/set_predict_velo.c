#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "physics.h"

// �\�����x�̌v�Z�̕⏕���s���֐�_1
void sub_set_predict_velo(int nNumClt, int nNumCltminusone);
// �\�����x�̌v�Z�̕⏕���s���֐�_2
void sub_sub_set_predict_velo(int nNumClt, int nNumCltminusone);

//BD///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// �\�����x�̌v�Z���s���֐�
void set_predict_velo(int nNumClt, int nNumClutOrigBranch)
{
	int i;

	// 0�Ԗڂ̍��̂̂Ƃ�
	if (nNumClt == 0)
	{
		for (i=0;i<6;++i)
		{
			clust[nNumClt].predict_velo[i] = 0.0;
		}
	}
	// �I�[�̍��̗̂\�����x�̌v�Z���s��
	else if(clust[nNumClt].terminal == TERMINAL)
	{
		sub_sub_set_predict_velo(nNumClt, nNumClutOrigBranch);
	}
	// ����̍��̗̂\�����x�̌v�Z���s��
	else if (clust[nNumClt-1].terminal == TERMINAL)
	{
		sub_set_predict_velo(nNumClt, nNumClutOrigBranch);
	}
	// �ʏ�̍��̗̂\�����x�̌v�Z���s��
	else
	{
		sub_set_predict_velo(nNumClt, nNumClt-1);
	}
}

// �\�����x�̌v�Z�̕⏕���s���֐�_1
void sub_set_predict_velo(int nNumClt, int nNumCltminusone)
{
	int i,j;

	// ���̗̂\�����x�̏��������s��
	for(i=0; i<6; ++i)
	{
		clust[nNumClt].predict_velo[i] = 0.0;
	}

	// �ϊ��s��̓]�u�s����悶��
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

// �\�������x�̌v�Z�̕⏕���s���֐�_2
void sub_sub_set_predict_velo(int nNumClt, int nNumCltminusone)
{
	int i,j;

	double TransMatrix[6][6];

	// ���̗̂\�����x�̏��������s��
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

	// �ϊ��s��̓]�u�s����悶��
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
