#include <stdio.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "force.h"
#include "MD.h"
#include "BD.h"

// BD /////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// ���݂̍\���ł̃^���p�N����
// spatial velocity �̌v�Z���s���֐�
void calc_sp_velo_cycle(int nNumClut)
{
	int i;

	// spatial velocity �̏��������s��
	for(i=0; i<6; ++i)
	{
		clust[nNumClut].sp_velo[i] = 0.0;
	}

	// spatial velocity �̌v�Z���s��
	sub_calc_sp_velo_cycle(nNumClut);
}

// ���݂̍\���ł̃^���p�N����
// spatial velocity �̌v�Z�̕⏕���s���֐�
void sub_calc_sp_velo_cycle(int nNumClut)
{
	int i;

	int alpha,alpha2;

	int nNumFreedom;

	double q[3];

	double OmegaOnThisCoord[3];
	double OmegaOnLaboCoord[3];

	// ���̂�spatial velocity �̏��������s��
	for(i=0; i<6; ++i)
	{
		clust[nNumClut].sp_velo[i] = 0.0;
	}


	// ���̂�2�ʊp�����x�̏��������s��
	clust[nNumClut].ddihedang[0] = 0.0;

////////////////////////////////////////////////////////////////////////
	// ���̂�2�ʊp���x�̌v�Z���s��_1
	for(i=0; i<6; ++i)
	{
//		clust[nNumClut].ddihedang[0]/*rad/s^2*/
//						-=   ABI[nNumClut].Kalman_Gain_Transpose[i]
//					    	*clust[nNumClut].predict_velo[i]/*rad/s^2*/;
	}

	// ���̂�2�ʊp���x�̌v�Z���s��_2
	//0410	clust[nNumClut].ddihedang[0]/*rad/s^2*/ 
	//				+=  zzz[nNumClut].nyu
	//				   /clust[nNumClut].friction_tensor_tra/*rad/s^2*/;
////////////////////////////////////////////////////////////////////////

	// spatial velocity �̌v�Z���s��_1
	for(i=0; i<6; ++i)
	{
		clust[nNumClut].sp_velo[i]/*rad/s^2*/
		                 +=   clust[nNumClut].predict_velo[i]/*rad/s^2*/;
	}

/*	// �N���X�^���W�n�ł̊p�����x��
	// ���������W�n�ł̊p�����x�̏�����
	for (alpha=0;alpha<3;++alpha)
	{
		OmegaOfThisCoord[alpha] = 0.0;
		OmegaOfAbsoCoord[alpha] = 0.0;
	}

	// ���̃N���X�^�̎��R�x�̎擾
	nNumFreedom = ABI[nNumClut].hingmat-1;

	// �N���X�^���W�n�ł̊p�����x�̎擾
	OmegaOfThisCoord[nNumFreedom] = clust[nNumClut].ddihedang[0];

	// ���������W�n�ł̊p�����x�̌v�Z
	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			OmegaOfAbsoCoord[alpha]
				+= clust[nNumClut].trans_A_to_CN[0][alpha2][alpha]
				  *OmegaOfThisCoord[alpha2];
		}
	}
*/

	// spatial acceleration �̌v�Z���s��_2
/*	for (alpha=0;alpha<3;++alpha)
	{
		clust[nNumClut].sp_velo[alpha]/*rad/s^2*/ 
/*		                 += OmegaOfAbsoCoord[alpha];/*rad/s^2*/
/*	}*/
 	// ���R�x�̑I��
	nNumFreedom = ABI[nNumClut].hingmat-1;

 	// �N���X�^���W�n�ł̊p���x�̏�����
	for (alpha=0;alpha<3;++alpha)
	{
		OmegaOnThisCoord[alpha] = 0.0;
		OmegaOnLaboCoord[alpha] = 0.0;
	}

 	// �N���X�^���W�n�ł̊p���x�̎擾
	OmegaOnThisCoord[nNumFreedom] 
	                 = clust[nNumClut].ddihedang[0];

 	// ���������W�n�ł̊p���x
/*	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			OmegaOnLaboCoord[alpha]
					 += clust[nNumClut].trans_A_to_CN[0][alpha2][alpha]
		               *OmegaOnThisCoord[alpha2];
		}
	}*/

	 // ���������W�n�ł̊p���x
	OmegaOnLaboCoord[0]
		=   clust[nNumClut].trans_A_to_CN[0][1][1]
		   *clust[nNumClut].trans_A_to_CN[0][2][0]
		   *clust[nNumClut].ddihedang[0]
		  - clust[nNumClut].trans_A_to_CN[0][1][0]
		   *clust[nNumClut].trans_A_to_CN[0][2][1]
		   *clust[nNumClut].ddihedang[0];

	OmegaOnLaboCoord[1]
		=   clust[nNumClut].trans_A_to_CN[0][0][1]
		   *clust[nNumClut].trans_A_to_CN[0][2][0]
		   *clust[nNumClut].ddihedang[0]
		  - clust[nNumClut].trans_A_to_CN[0][0][0]
		   *clust[nNumClut].trans_A_to_CN[0][2][1]
		   *clust[nNumClut].ddihedang[0];

	OmegaOnLaboCoord[2]
		=   clust[nNumClut].trans_A_to_CN[0][0][0]
		   *clust[nNumClut].trans_A_to_CN[0][1][1]
		   *clust[nNumClut].ddihedang[0]
		  - clust[nNumClut].trans_A_to_CN[0][0][1]
		   *clust[nNumClut].trans_A_to_CN[0][1][0]
		   *clust[nNumClut].ddihedang[0];

	// spatial velocity �̐ݒ���s��_2
	for (alpha=0;alpha<3;++alpha)
	{
		clust[nNumClut].sp_velo[alpha] += OmegaOnLaboCoord[alpha];
	}
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
