#include <stdio.h>
#include <math.h>
#include "gener.h"
#include "MD.h"
#include "ABA.h" // #include "ABA_multimer.h"

#include "EF.h"

// �Ǐ����W�n���������n�̕ϊ����s���֐�
void trans_CN_to_A(int nNumClut, int nNumClutOrigBranch)
{
	int i;
	int nNumAtomLoca;
	int nNumParent;

	nNumParent = clust[nNumClut].nNumClutOfParent-1;


	// 0�Ԗڂ̃N���X�^�ł͂Ȃɂ����Ȃ�
	if (nNumClut == 0)
	{
		;
	}
//	// �N���X�^���I�[�ł������Ƃ�
//	else if (clust[nNumClut].terminal == TERMINAL && nNumClut != prot.DOF-1)
//	{
//		// �Ǐ����W�n�̌��_�̌��q�ԍ��̎擾
//		nNumAtomLoca =   clust[nNumClutOrigBranch].num_atom_clust;
//
//		// �Ǐ����W�n���������n�̕ϊ�
//		sub_trans_CN_to_A(nNumClut,
//			              clust[nNumClut].num_atom_clust,
//			              nNumAtomLoca);
//	}
//	// �O�̃N���X�^���I�[�ł������Ƃ�
//	else if (clust[nNumClut-1].terminal == TERMINAL)
//	{
//		// �Ǐ����W�n�̌��_�̌��q�ԍ��̎擾
//		nNumAtomLoca =  clust[nNumClutOrigBranch].num_atom_clust
//		              + clust[nNumClut-1].num_atom_clust;
//
//		// �Ǐ����W�n���������n�̕ϊ�
//		sub_trans_CN_to_A(nNumClut,
//		             prot.num_atom-clust[nNumClut].origin_atom_a+1,
//		             nNumAtomLoca);
//	}	
//	// �ʏ�N���X�^�̂Ƃ�
//	else
//	{
//		// �Ǐ����W�n�̌��_�̌��q�ԍ��̎擾
//		nNumAtomLoca =  clust[nNumClut-1].num_atom_clust
//					   -clust[nNumClut-1].terminal_atom_a[0]
//					   +clust[nNumClut-1].origin_atom_a;
//
//		// �Ǐ����W�n���������n�̕ϊ�
//		sub_trans_CN_to_A(nNumClut,
//		              prot.num_atom-clust[nNumClut].origin_atom_a+1,
//		              nNumAtomLoca);
//	}
	else if (clust[nNumClut].terminal == TERMINAL && nNumClut != prot.DOF-1)
	{
		// �Ǐ����W�n�̌��_�̌��q�ԍ��̎擾
		nNumAtomLoca =   clust[nNumParent].num_atom_clust;

		// �Ǐ����W�n���������n�̕ϊ�
		sub_trans_CN_to_A(nNumClut,
			              clust[nNumClut].num_atom_clust,
			              nNumAtomLoca);
	}
	// �ʏ�N���X�^�̂Ƃ�
	else
	{
		// �Ǐ����W�n�̌��_�̌��q�ԍ��̎擾
		nNumAtomLoca =  clust[nNumParent].num_atom_clust
					   -clust[nNumParent].terminal_atom_a[0]
					   +clust[nNumParent].origin_atom_a;

		for (i=1;i<nNumClut-nNumParent;++i)
		{
			nNumAtomLoca += clust[nNumParent+i].num_atom_clust;
		}
//		nNumAtomLoca = clust[nNumClut].origin_atom_a;

		// �Ǐ����W�n���������n�̕ϊ�
		sub_trans_CN_to_A(nNumClut,
		              prot.num_atom-clust[nNumClut].origin_atom_a+1,
		              nNumAtomLoca);
	}

}

// �Ǐ����W�n���������n�̕ϊ����s���֐�
void sub_trans_CN_to_A(int nNumClt, int nNumAtom, int nNumAtomLoca)
{
	int i,j,k,alpha,alpha2;

	int n=0, angle;

	int nNumClut2;

	int n_delta_dihed;

	int nNumAtomAbsoOrig;

	double sn_delta_dihed;
	double cs_delta_dihed;

	double **Coord/*[MAXA][3]*/;
	double **Coord2/*[MAXA][3]*/;
	double Origin_Coord[3];

	double zaxis[3];
	double delta_dihed;

	double Rotation_Ele1;
	double Rotation_Ele2;
	double Rotation_Ele3;

	double Rotation_Mat1[3][3];
	double Rotation_Mat2[3][3];
	double Rotation_Mat3[3][3];

	double Rotation[3][3];

	FILE *outtest;

	Coord=(double **)gcemalloc(sizeof(double *)*prot.num_atom);
	Coord2=(double **)gcemalloc(sizeof(double *)*prot.num_atom);
	for (i=0;i<prot.num_atom;++i) {
	  Coord[i]=(double *)gcemalloc(sizeof(double)*3);
	  Coord2[i]=(double *)gcemalloc(sizeof(double)*3);
	}

	// ���X�e�b�v�ł̓�ʊp�̕ψ�
	delta_dihed = clust[nNumClt].now_deltadihedang[0];
//	delta_dihed = -0.1;
	//	delta_dihed = -clust[nNumClt].now_deltadihedang[0];
	if (nNumClt == 1 /*|| nNumClt ==*/ /*3*//*2*//*4*//*5*/)
	{
//		delta_dihed = -clust[nNumClt].now_deltadihedang[0];
	}

///////////////////////////////////////////////////////////
//	if (nNumClt == 3)
//	{
//		delta_dihed += 2.0*PI/10.0;
//	}
//	else
//	{
//		delta_dihed = 0.0;
//	}
///////////////////////////////////////////////////////////
/*	n_delta_dihed = (int)(delta_dihed/(2.0*PI));

	if (delta_dihed > 0.0)
	{
		delta_dihed -= (double)n_delta_dihed*2.0*PI;
	}
	else
	{
		delta_dihed += (double)n_delta_dihed*2.0*PI;
	}

	// ��ʊp�̕ψʂ̋K�i��_1
	if ( 0.0 >= delta_dihed )
	{
		for (;;)
		{
			delta_dihed += 2.0*PI;
			if ( delta_dihed >= 0.0 )
			{
				break;
			}
		}
	}

	// ��ʊp�̕ψʂ̋K�i��_2
	if ( delta_dihed >= 2.0*PI)
	{
		for (;;)
		{
			delta_dihed -= 2.0*PI;
			if ( delta_dihed <= 2.0*PI )
			{
				break;
			}
		}
	}
*/
	sn_delta_dihed = sin(delta_dihed*0.5);
	cs_delta_dihed = cos(delta_dihed*0.5);

	// ���_�Ǝn�_�̎擾
	nNumAtomAbsoOrig = clust[nNumClt].origin_atom_a-1;

	// ���_���W�̎擾
	for(alpha=0;alpha<3;++alpha)
	{
		Origin_Coord[alpha]=prot.coord[nNumAtomAbsoOrig][alpha];
	}

	// ��ʊp�̕ψʕ������Ǐ����W�̉�]
	for(i = 0; i < nNumAtom ; ++i)
	{
//		Coord[i][0] = clust[nNumClt].xoord_clust/*[0]*/[nNumAtomLoca+i][0];
//		Coord[i][1] = clust[nNumClt].xoord_clust/*[0]*/[nNumAtomLoca+i][1];
//		Coord[i][2] = clust[nNumClt].xoord_clust/*[0]*/[nNumAtomLoca+i][2];
		Coord[i][0] = prot.coord[nNumAtomAbsoOrig+i][0]-Origin_Coord[0];
		Coord[i][1] = prot.coord[nNumAtomAbsoOrig+i][1]-Origin_Coord[1];
		Coord[i][2] = prot.coord[nNumAtomAbsoOrig+i][2]-Origin_Coord[2];
	}

	for(alpha = 0; alpha < 3 ; ++alpha)
	{
		zaxis[alpha] = /*-*/clust[nNumClt].trans_A_to_CN[0][2][alpha];
	}

	for(i = 0; i < nNumAtom; ++i)
	{
		for(j = 0; j < 3 ; ++j)
		{
			Coord2[i][j] = 0.0;
		}
	}

	Rotation_Ele1 = cs_delta_dihed*cs_delta_dihed-sn_delta_dihed*sn_delta_dihed;
	Rotation_Ele2 = 2.0*sn_delta_dihed*sn_delta_dihed;
	Rotation_Ele3 = 2.0*cs_delta_dihed*sn_delta_dihed;

	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			Rotation_Mat1[alpha][alpha2] = ident[alpha][alpha2];
			Rotation_Mat2[alpha][alpha2] = zaxis[alpha]*zaxis[alpha2];
			Rotation_Mat3[alpha][alpha2] = 0.0;
		}
	}

	Rotation_Mat3[0][1] = -zaxis[2];
	Rotation_Mat3[0][2] =  zaxis[1];
	Rotation_Mat3[1][0] =  zaxis[2];
	Rotation_Mat3[1][2] = -zaxis[0];
	Rotation_Mat3[2][0] = -zaxis[1];
	Rotation_Mat3[2][1] =  zaxis[0];

	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			Rotation[alpha][alpha2] = Rotation_Ele1*Rotation_Mat1[alpha][alpha2]
									 +Rotation_Ele2*Rotation_Mat2[alpha][alpha2]
									 +Rotation_Ele3*Rotation_Mat3[alpha][alpha2];
		}
	}

	// �Ǐ����W�n���������n�̕ϊ�_�ϊ��s��̏�Z
	for(i = 0; i < nNumAtom; ++i)
	{
		for(j = 0; j < 3 ; ++j)
		{
			for(k = 0; k < 3 ; ++k)
			{
				Coord2[i][j] += Rotation[j][k]*Coord[i][k];
			}
		}
	}

	// �Ǐ����W�n���������n�̕ϊ�_���_�̈ړ�
	for(i = 0; i < nNumAtom ; ++i)
	{
		for(alpha=0 ;alpha<3; ++alpha)
		{
			prot.coord[nNumAtomAbsoOrig+i][alpha]
			                      = Coord2[i][alpha] + Origin_Coord[alpha];
		}
	}

	// ���̃X�e�b�v�̍��W�ł̓�ʊp�̐ݒ�
	clust[nNumClt].dihedang[0] += delta_dihed;

	// ��ʊp�̋K�i��_1
	if ( 0.0 >= clust[nNumClt].dihedang[0] )
	{
		for (;;)
		{
			clust[nNumClt].dihedang[0] += 2.0*PI;
			if ( clust[nNumClt].dihedang[0] >= 0.0 )
			{
				break;
			}
		}
	}

	// ��ʊp�̋K�i��_2
	if ( clust[nNumClt].dihedang[0] >= 2.0*PI)
	{
		for (;;)
		{
			clust[nNumClt].dihedang[0] -= 2.0*PI;
			if ( clust[nNumClt].dihedang[0] <= 2.0*PI )
			{
				break;
			}
		}
	}

	// ���̃X�e�b�v�̍��W�ł̋Ǐ����W�n�̌v�Z
	for(nNumClut2=nNumClt; nNumClut2<prot.DOF; ++nNumClut2)
	{
		trans_A_to_CN(nNumClut2);
	}

	if ((outtest=fopen("dihedang.out","a")) == NULL)
	{
		printf("in\n");
		exit(1);
	}

	fprintf(outtest, "%d %e \n",nNumClt, clust[nNumClt].dihedang[0]);

	fclose(outtest);
}
