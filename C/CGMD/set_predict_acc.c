#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "physics.h"

// �\�������x�̌v�Z�̕⏕���s���֐�_2
void sub_sub_set_predict_acc(int nNumClt, int nNumCltminusone);

// �\�������x�̌v�Z���s���֐�
void set_predict_acc(int nNumClt, int nNumClutOrigBranch)
{
	int i;
	int nNumClutOfParent;

	nNumClutOfParent = clust[nNumClt].nNumClutOfParent-1;

	// �N���X�^�����[�̂Ƃ�
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

// �\�������x�̌v�Z�̕⏕���s���֐�_1
void sub_set_predict_acc(int nNumClt, int nNumCltminusone)
{
	int i,j;

	// ���̗̂\�������x�̏��������s��
	for(i=0; i<6; ++i)
	{
		clust[nNumClt].predict_alpha[i] = 0.0;
	}

	// �ϊ��s��̓]�u�s����悶��
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

// �\�������x�̌v�Z�̕⏕���s���֐�_2
//void sub_sub_set_predict_acc(int nNumClt, int nNumCltminusone)
//{
//	int i,j;
//
//	double TransMatrix[6][6];
//
//	// ���̗̂\�������x�̏��������s��
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
//	// �ϊ��s��̓]�u�s����悶��
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
