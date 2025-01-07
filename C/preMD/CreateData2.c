#include <stdio.h>

#include "Vis_MD.h"

int TotalNumAtomType;

// �C���v�b�g�ɕK�v�ȃf�[�^�̍쐬���s���֐�
// ���W�f�[�^�̍쐬���s���֐�
int CreateDataCoord(elementoftable *element,
					int nNumAtomTotal);
// ���̃f�[�^�̍쐬���s���֐�
int CreateDataClust(elementoftable *element,
					int nNumAtomTotal,
					int nNumClutTotal);
// �͏���f�[�^�̍쐬���s���֐�
int CreateDataTop(elementoftable *element,
				  int nNumAtomTotal);
// ���q��̐����グ���s���֐�
int CountAtomID(elementoftable *element);

void CkeckNonBond(int nNumAtom);

double terminatecoord[3];

// �C���v�b�g�ɕK�v�ȃf�[�^�̍쐬�̐�����s���֐�
void CreateData(void)
{
	int nNumResidue;
	int AtomclustIDBEF;
	int HeadclustIDBEF=0;
s	int nNumAtom;
	int nNumAtom2;
	int nNumAtom3;

	int numTailatomonPro;
	int numTailatomonRes;
	int num_bond;
	int numTailclutonPro;
	int numTailclutonRes;
	int num_branch;

	char *nameresidue;
	elementoftable *elemnt;
	residuedata *dataofthisresidue;

	prot.nNumAtom = 0;
	prot.nNumClut = 0;

	for (nNumResidue=0;
		 nNumResidue<prot.nNumResidue;
		 ++nNumResidue)
	{
		// �c��f�[�^�̎擾���s��
		nameresidue = prot.Sequence[nNumResidue];
		elemnt = LURData(nameresidue,1,dataofthisresidue);

		// �c��Ԃ̌������s��
		// ���q
//		for (num_bond=0;;++num_bond)
//		{
//			if (prot.atm[numTailatomonPro].bond_interacting[num_bond]==0)
//			{
//				break;
//			}
//		}
		for (nNumAtom=0;
		 nNumAtom<elemnt->dataofthisresdiue.nNumAtom;
		 ++nNumAtom)
		{
			nNumAtom3 = prot.nNumAtom+nNumAtom;
			prot.atm[nNumAtom3].num_bond_interacting
			= elemnt->dataofthisresdiue.atm[numTailatomonRes].num_bond_interacting;
			for (nNumAtom2=0;
				nNumAtom2<prot.atm[nNumAtom3].num_bond_interacting;
				++nNumAtom2)
			{
				prot.atm[nNumAtom3].bond_interacting[nNumAtom2]
				= elemnt->dataofthisresdiue.atm[numTailatomonRes].bond_interacting[nNumAtom2];
			}
		}

		if (nNumResidue == 0)
		{
			numTailatomonPro = elemnt->dataofthisresdiue.TailAtom;
		}
		numTailatomonRes = elemnt->dataofthisresdiue.TailAtom;
		num_bond = elemnt->dataofthisresdiue.atm[numTailatomonRes].num_bond_interacting;
		prot.atm[numTailatomonPro].bond_interacting[num_bond]
		= elemnt->dataofthisresdiue.HeadAtom + prot.nNumAtom;
		prot.atm[numTailatomonRes].num_bond_interacting += 1;

		// �c��Ԃ̌������s��
		// ����
		if (nNumResidue == 0)
		{
			numTailclutonPro = elemnt->dataofthisresdiue.TailClut + prot.nNumClut;
		}
		num_branch = prot.clt[numTailclutonPro].nNumBranch;
		prot.clt[numTailclutonPro].nNumChildClust[num_branch+1]
		= elemnt->dataofthisresdiue.HeadClut + prot.nNumClut;
		prot.clt[numTailclutonPro].nNumBranch += 1;

		// ���q��𐔂��グ
		//		TotalNumAtomType = CountAtomID();

		// ���W�f�[�^�̍쐬���s��
//		CreateDataCoord(elemnt, prot.nNumAtom);
		// ���̃f�[�^�̍쐬���s��
		CreateDataClust(elemnt, prot.nNumAtom, prot.nNumClut);
		// �͏���f�[�^�̍쐬���s��
//		CreateDataTop(elemnt, prot.nNumAtom);

		// �^���p�N�����̑����q���̌v�Z
		prot.nNumAtom += elemnt->dataofthisresdiue.nNumAtom;
		// �^���p�N�����̑����̐��̌v�Z
		if (nNumResidue < prot.nNumResidue-1)
		{
			prot.nNumClut += elemnt->dataofthisresdiue.nNumClut-1;
		}
		else
		{
			prot.nNumClut += elemnt->dataofthisresdiue.nNumClut;
		}
		numTailatomonPro = elemnt->dataofthisresdiue.TailAtom + prot.nNumAtom;
		numTailclutonPro = elemnt->dataofthisresdiue.TailClut + prot.nNumClut;
		numTailclutonRes = elemnt->dataofthisresdiue.TailClut;
		//		num_branch = elemnt->dataofthisresdiue.atm[numTailclutonRes].nNum_Branch;

	}

	prot.nNumAtomType = TotalNumAtomType;
	prot.clt[prot.nNumClut-1].nNumBranch = 0;

//	for (nNumAtom=0;
//		 nNumAtom<prot.nNumAtom;
//		 ++nNumAtom)
//	{
//		for (nNumAtom2=0;
//			 /*nNumAtom2<prot.atm[nNumAtom]*/;
//		 	++nNumAtom2)
//		{
//			// �񌋍������ݍ�p���Ȃ��g
////			prot.atm[nNumAtom].not_interacting[];
//		}
//	}


}

// ���W�f�[�^�̍쐬���s���֐�
int CreateDataCoord(elementoftable *element, int nNumAtomTotal)
{
	int nNumAtom;
	int alpha;

	if (nNumAtomTotal == 0)
	{
		for (alpha=0;alpha<3;++alpha)
		{
			terminatecoord[alpha] = 0.0;
		}
	}

	for (nNumAtom = 0;
		 nNumAtom < element->dataofthisresdiue.nNumAtom-1;
		 ++nNumAtom)
	{
		for (alpha=0;alpha<3;++alpha)
		{
			// ���̌��q�̍��W�𑫂�
			prot.atm[nNumAtom+nNumAtomTotal].coord[alpha]
			=   element->dataofthisresdiue.atm[nNumAtom].coord[alpha]
			  + terminatecoord[alpha];
		}
	}

	for (alpha=0;alpha<3;++alpha)
	{
		terminatecoord[alpha]
			+= element->dataofthisresdiue.atm[nNumAtom].coord[alpha];
	}

//	void TransLocalToAbso (int nNumAtom1, int nNumAtom2,
//                         int nSumNumAtom, double delta_dihed)

}

// ���̃f�[�^�̍쐬���s���֐�
int CreateDataClust(elementoftable *element,
					int nNumAtomTotal,
					int nNumClutTotal)
{

	int nNumClust;
	int nNumBranch;
	int numChildclust;

	// ���̌��_�̌��q�ԍ�
	if (nNumClutTotal == 0)
	{
		prot.clt[0].origin = 1;
		prot.clt[0].nNumBranch
		= element->dataofthisresdiue.clt[0].nNumBranch;
		// ���̂��I�_���ǂ���
		prot.clt[0].termflag
		= element->dataofthisresdiue.clt[0].termflag;
		// ���̂̌��q��
		prot.clt[0].nNumAtom
		= element->dataofthisresdiue.clt[0].nNumAtom;
		// ���̂̎��R�x
		prot.clt[0].nDegOfFoo
		= 3;
		// ���̏I�_�̌��q�ԍ�
		prot.clt[nNumClust+nNumClutTotal].term[nNumBranch]
		= 1;
	}

	for (nNumClust = 0;
		 nNumClust < element->dataofthisresdiue.nNumClut;
		 ++nNumClust)
	{
		if (element->dataofthisresdiue.clt[nNumClust].origin != -1)
		{
			// ���̌��_�̌��q�ԍ�
			prot.clt[nNumClust+nNumClutTotal].origin
				= element->dataofthisresdiue.clt[nNumClust].origin + nNumAtomTotal;
			// ���̂̎}�̐�
			prot.clt[nNumClust+nNumClutTotal].nNumBranch
				= element->dataofthisresdiue.clt[nNumClust].nNumBranch;
			// ���̂��I�_���ǂ���
			prot.clt[nNumClust+nNumClutTotal].termflag
				= element->dataofthisresdiue.clt[nNumClust].termflag;
			// ���̂̌��q��
			prot.clt[nNumClust+nNumClutTotal].nNumAtom
				= element->dataofthisresdiue.clt[nNumClust].nNumAtom;
			// ���̂̎��R�x
			prot.clt[nNumClust+nNumClutTotal].nDegOfFoo
				= 3;

			for (nNumBranch=0;
				 nNumBranch<prot.clt[nNumClust].nNumBranch;
					++nNumBranch)
			{
				// ���̏I�_�̌��q�ԍ�
				prot.clt[nNumClust+nNumClutTotal].term[nNumBranch]
				=   element->dataofthisresdiue.clt[nNumClust].term[nNumBranch]
				   + nNumAtomTotal;
			}
		}
	}
	for (nNumClust = 0;
		 nNumClust < element->dataofthisresdiue.nNumClut;
		 ++nNumClust)
	{
		for (nNumBranch=0;
			 nNumBranch<prot.clt[nNumClust].nNumBranch;
				++nNumBranch)
		{
			// "�e"�� OK �Ȃ̂ŁA"�q"���쐬
//			numChildclust = elemnt->dataofthisresdiue->clt[nNumClut+nNumClutTotal].nNumChildClut;
//			elemnt->dataofthisresdiue->clt[nNumClut+nNumClutTotal].nNumChildClut[nNumBranch] = numChildclust;
		}
	}

}

// �͏���f�[�^�̍쐬���s���֐�
int CreateDataTop(elementoftable *element, int nNumAtomTotal)
{
	int nNumAtom;
	int nNumClust;
	int nNumDihed;
	int nNumAtomoforigin;
	int nNumClutTotal;

	// ���q��𐔂��グ
	TotalNumAtomType = CountAtomID();

	// ��ʊp�f�[�^�̍\�z//////////////////////////////////////////////////////
	for (nNumClust = 0;
		 nNumClust < element->dataofthisresdiue.nNumClut;
		 ++nNumClust)
	{
//		for (nNumBondInt;;++nNumBondInt)
//		{
//			if ((nNUmAtomOfBond
//				 = element->dataofthisresdiue.atm[nNumAtom].bond_interacting[nNumBondInt])
//				!= NULL)
//			{
//				prot.atm[nNumAtomTotal+nNumAtom].bond_interacting[nNumBondInt]
//				= nNUmAtomOfBond
//				 +nNumAtomTotal;
//			}
//			else
//			{
//				break;
//			}
//		}
	}

	nNumAtomoforigin = prot.clt[nNumClutTotal].origin;
//	prot.atm[nNumAtomoforigin].bond_interacting[nNumBondInt] = nNumAtomTotal;

//	nNumAtomofterm = prot.clt[nNumClust+nNumClutTotal].term[0];
//	prot.atm[nNumAtomofterm].bond_interacting[nNumBondInt] = nNumAtom+nNumAtomTotal;

//	for (nNumDihed=0;;++nNumDihed)
//	{
//		prot.dihedang[nNumDihed].atomnum[1] = nNumAtomoforigin;
//		prot.dihedang[nNumDihed].atomnum[2] = nNumAtomofterm;
//		prot.dihedang[nNumDihed].atomnum[0] = prot.atm[nNumAtomofterm].bond_interacting[nNumBondInt];
//		prot.dihedang[nNumDihed].atomnum[3] = prot.atm[nNumAtomofterm].bond_interacting[nNumBondInt];
//	}
	///////////////////////////////////////////////////////////////////////////

	// ���q�f�[�^�̍\�z//////////////////////////////////////////////////////
	for (nNumAtom=0;
		 nNumAtom<element->dataofthisresdiue.nNumAtom;
		 ++nNumAtom)
	{
//		// gamma
//		prot.atm[nNumAtomTotal+nNumAtom].gamma
//		= element->dataofthisresdiue.atm[nNumAtom].atomID.gamma;
//		// eata
//		prot.atm[nNumAtomTotal+nNumAtom].eata
//		= element->dataofthisresdiue.atm[nNumAtom].atomID.eata;
	}
	///////////////////////////////////////////////////////////////////////////

	for (nNumAtom=0;
		 nNumAtom<element->dataofthisresdiue.nNumAtom;
		 ++nNumAtom)
	{
		CkeckNonBond(nNumAtom);
	  //		for (nNumAtom2=0;
	  //			 ;nNumAtom2 < element->dataofthisresdiue.atm[nNumAtom].num_bond;
	  //			 ;++nNumAtom2)
	  //		{
	  //			prot.atm[nNumAtomTotal+nNumAtom].not_interacting[nNumAtom2]
	  //		= element->dataofthisresdiue.atm[nNumAtom].bond_interacting[nNumAtom2];
	  //		prot.atm[nNumAtomTotal+nNumAtom].o_f_interacting[nNumAtom2]
	  //		= element->dataofthisresdiue.atm[nNumAtom];
	  //	}
	}


}

// ���q��̐����グ���s���֐�
int CountAtomID(void)
{
	int TotalNumAtomType=0;
	int ThisAtomType;
	int numAtom;
	int flag = 0;
	int ON = 1;
	int OFF = 0;
	int list[100/*MAXNAT*/];
	int i,j;

	for (numAtom=0;numAtom<prot.nNumAtom;++numAtom)
	{
		ThisAtomType = prot.atm[numAtom].atomID.atomnum;
		for (j=0;list[j]!=NULL;++j)
		{
			if (ThisAtomType != list[j])
			{
				flag = ON;
			}
		}
		if (flag == ON)
		{
			list[j] = ThisAtomType;
			++TotalNumAtomType;
		}
		flag = OFF;
	}

	return TotalNumAtomType;
}

void CkeckNonBond(int nNumAtom)
{
	int i;
	int num=0;
	int atomnum[3];

	prot.atm[nNumAtom].not_interacting[num]=nNumAtom;
	// n+1 �Ԗڂ̌��q
	for (i=0;i<prot.atm[nNumAtom].num_bond_interacting;++i)
	{
		++num;
		prot.atm[nNumAtom].not_interacting[num]
		=prot.atm[nNumAtom].bond_interacting[i];
		atomnum[0]=prot.atm[nNumAtom].bond_interacting[i];
		// n+2 �Ԗڂ̌��q
		for (i=0;i<prot.atm[atomnum[0]].num_bond_interacting;++i)
		{
			++num;
			prot.atm[nNumAtom].not_interacting[num]
			=prot.atm[nNumAtom].bond_interacting[i];
			atomnum[1]=prot.atm[nNumAtom].bond_interacting[i];
			// n+3 �Ԗڂ̌��q
			for (i=0;i<prot.atm[atomnum[1]].num_bond_interacting;++i)
			{
				++num;
				prot.atm[nNumAtom].not_interacting[num]
				=prot.atm[nNumAtom].bond_interacting[i];
				atomnum[2]=prot.atm[nNumAtom].bond_interacting[i];
				// n+4 �Ԗڂ̌��q
				for (i=0;i<prot.atm[atomnum[2]].num_bond_interacting;++i)
				{
					++num;
					prot.atm[nNumAtom].not_interacting[num]
					=prot.atm[nNumAtom].bond_interacting[i];
//					prot.atm[nNumAtom].o_f_interacting[num]
//					=prot.atm[nNumAtom].bond_interacting[i];
				}
			}
		}
	}
}
