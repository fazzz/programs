#include <stdio.h>
#include <math.h>
#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"

// �������n���Ǐ����W�n�̕ϊ����s���֐�
void trans_A_to_CN(int nNumClut)
{
	int nNumClutOrigBranch;
	int nNumClutOfParent;
	int nNumClut2;
	int num_atom;
	int num,nNumClutdummy;
	int flag;

	nNumClutOfParent = clust[nNumClut].nNumClutOfParent-1;

	if (clust[nNumClut].num_branch > 1)
	{
		nNumClutOrigBranch = nNumClut;
	}

	// 0�Ԗڂ̃N���X�^�ł́A���_�̈ړ����s��
	if (nNumClut == 0)
	{
		sub_trans_A_to_CN(0, -1,0, 
					       prot.num_atom);
		sub_trans_A_to_CN_Initial();
	}
//	// �N���X�^���I�[�ł������Ƃ�
//	else if(clust[nNumClut].terminal == TERMINAL && nNumClut != prot.DOF-1)
//	{
//		sub_trans_A_to_CN(nNumClut, nNumClutOrigBranch,0, 
//					       clust[nNumClutOrigBranch].num_atom_clust
//					      +clust[nNumClut].num_atom_clust);
//	}
//	// �O�̃N���X�^���I�[�ł������Ƃ�
//	else if(clust[nNumClut-1].terminal == TERMINAL)
//	{
//		sub_trans_A_to_CN(nNumClut, nNumClutOrigBranch, 0, 
//						   prot.num_atom-clust[nNumClutOrigBranch].terminal_atom_a[0]+1);
//	}
//	// �ʏ�N���X�^�̂Ƃ�
//	else
//	{
//		sub_trans_A_to_CN(nNumClut, nNumClut-1, 0,
//		                   prot.num_atom-clust[nNumClut-1].terminal_atom_a[0]+1);
//	}
	// �N���X�^���I�[�ł������Ƃ�
	else if(clust[nNumClut].terminal == TERMINAL && nNumClut != prot.DOF-1)
	{
	  num=0;
	  //	  sub_trans_A_to_CN(nNumClut, nNumClutOfParent,0, 
	  //		    clust[nNumClutOfParent].num_atom_clust
	  //		    +clust[nNumClut].num_atom_clust);
	  for (nNumClutdummy=nNumClutOfParent;nNumClutdummy<=nNumClut;++nNumClutdummy) {
	    num+=clust[nNumClutdummy].num_atom_clust;
	  }
	  sub_trans_A_to_CN(nNumClut, nNumClutOfParent,0, 
			    num);
	}
	else if(clust[nNumClut].terminal == BRANCH)
	{
	        num_atom = clust[nNumClutOfParent].num_atom_clust;flag=1;
	        for (nNumClut2=nNumClut;/*clust[nNumClut2].terminal != TERMINAL*/
		     (clust[nNumClut2].terminal != TERMINAL && flag ==1) ||  (clust[nNumClut2].terminal == TERMINAL && clust[nNumClut2].nNumClutOfParent-1 == nNumClut);
		     ++nNumClut2)
		{
		  if (clust[nNumClut2].terminal == TERMINAL)
		    flag = 0;
		  num_atom+=clust[nNumClut2].num_atom_clust;
		}
		//		num_atom+=clust[nNumClut2].num_atom_clust;
		sub_trans_A_to_CN(nNumClut, nNumClutOfParent,0, num_atom);
	}
	// �O�̃N���X�^���I�[�ł������Ƃ�
//	else if(clust[nNumClut-1].terminal == TERMINAL)
//	{
//		sub_trans_A_to_CN(nNumClut, nNumClutOrigBranch, 0, 
//						   prot.num_atom-clust[nNumClutOrigBranch].terminal_atom_a[0]+1);
//	}
	// �ʏ�N���X�^�̂Ƃ�
	else
	{
		sub_trans_A_to_CN(nNumClut, nNumClutOfParent, 0,
		                   prot.num_atom-clust[nNumClutOfParent].terminal_atom_a[0]+1);
	}
}

// �������n���Ǐ����W�n�̕ϊ��̕⏕���s���֐�
void sub_trans_A_to_CN(int nNumCltTar, int nNumCltCoo,
                       int nNumBod, int nNumAtom)
{
	int i;

	int alpha;

	int nNmAtomOfCN_A;
	int nNmAtomOfCN_1_A;

	double CN_A[3];
	double HN_A[3];
	double CN_1_A[3];
	double mat[MAXA][3];

	double ii[3];
	double jj[3];
	double kk[3];

	double jj_x;
	double jj_y;
	double jj_z;

	double DisOfCN_1_CN=0.0;
	double DisOfHN_CN=0.0;
	double SnCN_1_CN_HN=0.0;
	double CsCN_1_CN_HN=0.0;

	clust[nNumCltTar].num_xoord_a = nNumAtom;

	// ���q�ԍ��̎擾
	nNmAtomOfCN_A = clust[nNumCltTar].origin_atom_a-1;
	if (nNumCltCoo != -1)
	{
		nNmAtomOfCN_1_A = clust[nNumCltCoo].terminal_atom_a[0]-1;
	}
	else
	{
		nNmAtomOfCN_1_A = 0;
	}

	clust[nNumCltTar].origin_xoord_a = nNmAtomOfCN_A-nNmAtomOfCN_1_A;

	// �������n�ł̍��W�̎擾
	for(alpha=0;alpha<3;++alpha)
	{
		CN_A[alpha]/*A*/=prot.coord[nNmAtomOfCN_A][alpha]/*A*/;
		HN_A[alpha]/*A*/=prot.coord[nNmAtomOfCN_A+1][alpha]/*A*/;
		if (nNumCltCoo != -1)
		{
			CN_1_A[alpha]/*A*/=prot.coord[nNmAtomOfCN_1_A][alpha]/*A*/;
		}
		else
		{
			CN_1_A[alpha]/*A*/=0.0/*A*/;
		}
	}

	for(alpha=0;alpha<3;++alpha)
	{
		// ^0k1�̌v�Z_1
		kk[alpha]/*A*/ = CN_1_A[alpha]-CN_A[alpha]/*A*/;
		// |CN_1-CN|�̌v�Z_1
		DisOfCN_1_CN
		 += (CN_1_A[alpha]-CN_A[alpha])*(CN_1_A[alpha]-CN_A[alpha]);
		// |CN_1-CN|�̌v�Z_1
		DisOfHN_CN
		 += (HN_A[alpha]-CN_A[alpha])*(HN_A[alpha]-CN_A[alpha]);
		// �p(CN_1-CN-HN)�̌v�Z_1
		CsCN_1_CN_HN += (CN_1_A[alpha]-CN_A[alpha])*(HN_A[alpha]-CN_A[alpha]);
	}

	// |CN_1-CN|�̌v�Z_2
	DisOfCN_1_CN/*A*/ = sqrt(DisOfCN_1_CN)/*A*/;
	// |CN_1-CN|�̌v�Z_2
	DisOfHN_CN/*A*/ = sqrt(DisOfHN_CN)/*A*/;
	// �p(CN_1-CN-HN)�̌v�Z_2
	CsCN_1_CN_HN = /*-*/CsCN_1_CN_HN/(DisOfCN_1_CN*DisOfHN_CN);
	SnCN_1_CN_HN = 1.0-CsCN_1_CN_HN*CsCN_1_CN_HN;
	SnCN_1_CN_HN = sqrt(SnCN_1_CN_HN);

	// ^0k1�̌v�Z_2
	for(alpha=0;alpha<3;++alpha)
	{
		kk[alpha]=kk[alpha]/DisOfCN_1_CN;
	}

	// ^0j1�̌v�Z
	jj[0]=((CN_1_A[1]-CN_A[1])*(HN_A[2]-CN_A[2])
		  -(CN_1_A[2]-CN_A[2])*(HN_A[1]-CN_A[1]))
		  /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);
	jj[1]=((CN_1_A[2]-CN_A[2])*(HN_A[0]-CN_A[0])
	      -(CN_1_A[0]-CN_A[0])*(HN_A[2]-CN_A[2]))
	      /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);
	jj[2]=((CN_1_A[0]-CN_A[0])*(HN_A[1]-CN_A[1])
	      -(CN_1_A[1]-CN_A[1])*(HN_A[0]-CN_A[0]))
	      /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);

//	jj[0]= -jj[0];
//	jj[1]= -jj[1];
//	jj[2]= -jj[2];

//	// ^0j1�̌v�Z
//	jj[0]=-( (CN_1_A[1]-CN_A[1])*(HN_A[2]-CN_A[2])
//		    -(CN_1_A[2]-CN_A[2])*(HN_A[1]-CN_A[1]))
//		  /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);
//	jj[1]=-( (CN_1_A[2]-CN_A[2])*(HN_A[0]-CN_A[0])
//	        -(CN_1_A[0]-CN_A[0])*(HN_A[2]-CN_A[2]))
//	      /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);
//	jj[2]=-( (CN_1_A[0]-CN_A[0])*(HN_A[1]-CN_A[1])
//	        -(CN_1_A[1]-CN_A[1])*(HN_A[0]-CN_A[0]))
//	      /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);

	// ^0i1�̌v�Z
//	ii[0]=(jj[1]*(CN_1_A[2]-CN_A[2])-jj[2]*(CN_1_A[1]-CN_A[1]))/DisOfCN_1_CN;
//	ii[1]=(jj[2]*(CN_1_A[0]-CN_A[0])-jj[0]*(CN_1_A[2]-CN_A[2]))/DisOfCN_1_CN;
//	ii[2]=(jj[0]*(CN_1_A[1]-CN_A[1])-jj[1]*(CN_1_A[0]-CN_A[0]))/DisOfCN_1_CN;

	ii[0]=jj[1]*kk[2]-jj[2]*kk[1];
	ii[1]=jj[2]*kk[0]-jj[0]*kk[2];
	ii[2]=jj[0]*kk[1]-jj[1]*kk[0];

//	// ^0i1�̌v�Z
//	ii[0]=-(jj[1]*(CN_1_A[2]-CN_A[2])-jj[2]*(CN_1_A[1]-CN_A[1]))/DisOfCN_1_CN;
//	ii[1]=-(jj[2]*(CN_1_A[0]-CN_A[0])-jj[0]*(CN_1_A[2]-CN_A[2]))/DisOfCN_1_CN;
//	ii[2]=-(jj[0]*(CN_1_A[1]-CN_A[1])-jj[1]*(CN_1_A[0]-CN_A[0]))/DisOfCN_1_CN;

//	// ^0j1�̌v�Z
//	ii[0]=((CN_1_A[1]-CN_A[1])*(HN_A[2]-CN_A[2])
//		  -(CN_1_A[2]-CN_A[2])*(HN_A[1]-CN_A[1]))
//		  /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);
//	ii[1]=((CN_1_A[2]-CN_A[2])*(HN_A[0]-CN_A[0])
//	      -(CN_1_A[0]-CN_A[0])*(HN_A[2]-CN_A[2]))
//	      /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);
//	ii[2]=((CN_1_A[0]-CN_A[0])*(HN_A[1]-CN_A[1])
//	      -(CN_1_A[1]-CN_A[1])*(HN_A[0]-CN_A[0]))
//	      /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);
//
////	// ^0j1�̌v�Z
////	jj[0]=-( (CN_1_A[1]-CN_A[1])*(HN_A[2]-CN_A[2])
////		    -(CN_1_A[2]-CN_A[2])*(HN_A[1]-CN_A[1]))
////		  /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);
////	jj[1]=-( (CN_1_A[2]-CN_A[2])*(HN_A[0]-CN_A[0])
////	        -(CN_1_A[0]-CN_A[0])*(HN_A[2]-CN_A[2]))
////	      /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);
////	jj[2]=-( (CN_1_A[0]-CN_A[0])*(HN_A[1]-CN_A[1])
////	        -(CN_1_A[1]-CN_A[1])*(HN_A[0]-CN_A[0]))
////	      /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);
//
//	// ^0i1�̌v�Z
//	jj[0]=(ii[1]*(CN_1_A[2]-CN_A[2])-ii[2]*(CN_1_A[1]-CN_A[1]))/DisOfCN_1_CN;
//	jj[1]=(ii[2]*(CN_1_A[0]-CN_A[0])-ii[0]*(CN_1_A[2]-CN_A[2]))/DisOfCN_1_CN;
//	jj[2]=(ii[0]*(CN_1_A[1]-CN_A[1])-ii[1]*(CN_1_A[0]-CN_A[0]))/DisOfCN_1_CN;
//
////	// ^0i1�̌v�Z
////	ii[0]=-(jj[1]*(CN_1_A[2]-CN_A[2])-jj[2]*(CN_1_A[1]-CN_A[1]))/DisOfCN_1_CN;
////	ii[1]=-(jj[2]*(CN_1_A[0]-CN_A[0])-jj[0]*(CN_1_A[2]-CN_A[2]))/DisOfCN_1_CN;
////	ii[2]=-(jj[0]*(CN_1_A[1]-CN_A[1])-jj[1]*(CN_1_A[0]-CN_A[0]))/DisOfCN_1_CN;

	// ^0j1�̌v�Z
//	ii[0]=((CN_1_A[1]-CN_A[1])*(HN_A[2]-CN_A[2])
//		  -(CN_1_A[2]-CN_A[2])*(HN_A[1]-CN_A[1]))
//		  /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);
//	ii[1]=((CN_1_A[2]-CN_A[2])*(HN_A[0]-CN_A[0])
//	      -(CN_1_A[0]-CN_A[0])*(HN_A[2]-CN_A[2]))
//	      /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);
//	ii[2]=((CN_1_A[0]-CN_A[0])*(HN_A[1]-CN_A[1])
//	      -(CN_1_A[1]-CN_A[1])*(HN_A[0]-CN_A[0]))
//	      /(DisOfCN_1_CN*DisOfHN_CN*SnCN_1_CN_HN);

	// ^0i1�̌v�Z
//	jj[0]=(ii[1]*(CN_1_A[2]-CN_A[2])-ii[2]*(CN_1_A[1]-CN_A[1]))/DisOfCN_1_CN;
//	jj[1]=(ii[2]*(CN_1_A[0]-CN_A[0])-ii[0]*(CN_1_A[2]-CN_A[2]))/DisOfCN_1_CN;
//	jj[2]=(ii[0]*(CN_1_A[1]-CN_A[1])-ii[1]*(CN_1_A[0]-CN_A[0]))/DisOfCN_1_CN;

	// ���W�n�̕ϊ�_���_�̏d�ˍ��킹
	for(i=0; i<nNumAtom; ++i)
	{
		for(alpha=0;alpha<3;++alpha)
		{
			mat[i][alpha] = prot.coord[nNmAtomOfCN_1_A+i][alpha]/*A*/
			              -CN_A[alpha]/*A*/;
		}
	}

	// ���W�n�̕ϊ��s��̍쐬
	for(alpha=0;alpha<3;++alpha)
	{
		clust[nNumCltTar].trans_A_to_CN[nNumBod][0][alpha]=ii[alpha];
		clust[nNumCltTar].trans_A_to_CN[nNumBod][1][alpha]=jj[alpha];
		clust[nNumCltTar].trans_A_to_CN[nNumBod][2][alpha]=kk[alpha];
	}

//	// ���W�n�̕ϊ��s��̍쐬
//	for(alpha=0;alpha<3;++alpha)
//	{
//		clust[nNumCltTar].trans_A_to_CN[nNumBod][alpha][0]=ii[alpha];
//		clust[nNumCltTar].trans_A_to_CN[nNumBod][alpha][1]=jj[alpha];
//		clust[nNumCltTar].trans_A_to_CN[nNumBod][alpha][2]=kk[alpha];
//	}

	// ���W�n�̕ϊ�_�ϊ��s��̏�Z
	for(i=0; i < nNumAtom; ++i)
	{
		clust[nNumCltTar].xoord_clust[nNumBod][i][0]
		 =   clust[nNumCltTar].trans_A_to_CN[nNumBod][0][0]*mat[i][0]
		   + clust[nNumCltTar].trans_A_to_CN[nNumBod][0][1]*mat[i][1]
		   + clust[nNumCltTar].trans_A_to_CN[nNumBod][0][2]*mat[i][2];

		clust[nNumCltTar].xoord_clust[nNumBod][i][1]
		 =   clust[nNumCltTar].trans_A_to_CN[nNumBod][1][0]*mat[i][0]
		   + clust[nNumCltTar].trans_A_to_CN[nNumBod][1][1]*mat[i][1]
		   + clust[nNumCltTar].trans_A_to_CN[nNumBod][1][2]*mat[i][2];

		clust[nNumCltTar].xoord_clust[nNumBod][i][2]
		 =   clust[nNumCltTar].trans_A_to_CN[nNumBod][2][0]*mat[i][0]
		   + clust[nNumCltTar].trans_A_to_CN[nNumBod][2][1]*mat[i][1]
		   + clust[nNumCltTar].trans_A_to_CN[nNumBod][2][2]*mat[i][2];
	}
}

// 0�Ԗڂ̍��̂̎������n���Ǐ����W�n�̕ϊ��̕⏕���s���֐�
void sub_trans_A_to_CN_Initial(void)
{
	int i;
	int alpha, alpha2;

	// �P�ʍs��̐ݒ���s��
	ident[0][0]=1.0;
	ident[0][1]=0.0;
	ident[0][2]=0.0;
	ident[1][0]=0.0;
	ident[1][1]=1.0;
	ident[1][2]=0.0;
	ident[2][0]=0.0;
	ident[2][1]=0.0;
	ident[2][2]=1.0;
	ident[0][3]=0.0;
	ident[0][4]=0.0;
	ident[0][5]=0.0;
	ident[1][3]=0.0;
	ident[1][4]=0.0;
	ident[1][5]=0.0;
	ident[2][3]=0.0;
	ident[2][4]=0.0;
	ident[2][5]=0.0;
	ident[3][0]=0.0;
	ident[3][1]=0.0;
	ident[3][2]=0.0;
	ident[3][3]=1.0;
	ident[3][4]=0.0;
	ident[3][5]=0.0;
	ident[4][0]=0.0;
	ident[4][1]=0.0;
	ident[4][2]=0.0;
	ident[4][3]=0.0;
	ident[4][4]=1.0;
	ident[4][5]=0.0;
	ident[5][0]=0.0;
	ident[5][1]=0.0;
	ident[5][2]=0.0;
	ident[5][3]=0.0;
	ident[5][4]=0.0;
	ident[5][5]=1.0;

	// ���W�n�̕ϊ�
	for (i=0;i<prot.num_atom;++i)
	{
		for (alpha=0;alpha<3;++alpha)
		{
//			prot.coord[i][alpha] = clust[0].xoord_clust/*[0]*/[i][alpha];
		}
	}

	for(alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			clust[0].trans_A_to_CN[0][alpha][alpha2]=0.0;
		}
	}

	// ���W�n�̕ϊ��s��̍쐬
	for(alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			clust[0].trans_A_to_CN[0][alpha][alpha2]=ident[alpha][alpha2];
		}
	}

	clust[0].trans_A_to_CN[0][0][0]=1.0;
	clust[0].trans_A_to_CN[0][0][1]=0.0;
	clust[0].trans_A_to_CN[0][0][2]=0.0;
	clust[0].trans_A_to_CN[0][1][0]=0.0;
	clust[0].trans_A_to_CN[0][1][1]=1.0;
	clust[0].trans_A_to_CN[0][1][2]=0.0;
	clust[0].trans_A_to_CN[0][2][0]=0.0;
	clust[0].trans_A_to_CN[0][2][1]=0.0;
	clust[0].trans_A_to_CN[0][2][2]=1.0;
}
