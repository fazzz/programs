#include <stdio.h>

#include "Vis_MD.h"

int TotalNumAtomType;

// インプットに必要なデータの作成を行う関数
// 座標データの作成を行う関数
int CreateDataCoord(elementoftable *element,
					int nNumAtomTotal);
// 剛体データの作成を行う関数
int CreateDataClust(elementoftable *element,
					int nNumAtomTotal,
					int nNumClutTotal);
// 力場情報データの作成を行う関数
int CreateDataTop(elementoftable *element,
				  int nNumAtomTotal);
// 原子種の数え上げを行う関数
int CountAtomID(elementoftable *element);

void CkeckNonBond(int nNumAtom);

double terminatecoord[3];

// インプットに必要なデータの作成の制御を行う関数
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
		// 残基データの取得を行う
		nameresidue = prot.Sequence[nNumResidue];
		elemnt = LURData(nameresidue,1,dataofthisresidue);

		// 残基間の結合を行う
		// 原子
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

		// 残基間の結合を行う
		// 剛体
		if (nNumResidue == 0)
		{
			numTailclutonPro = elemnt->dataofthisresdiue.TailClut + prot.nNumClut;
		}
		num_branch = prot.clt[numTailclutonPro].nNumBranch;
		prot.clt[numTailclutonPro].nNumChildClust[num_branch+1]
		= elemnt->dataofthisresdiue.HeadClut + prot.nNumClut;
		prot.clt[numTailclutonPro].nNumBranch += 1;

		// 原子種を数え上げ
		//		TotalNumAtomType = CountAtomID();

		// 座標データの作成を行う
//		CreateDataCoord(elemnt, prot.nNumAtom);
		// 剛体データの作成を行う
		CreateDataClust(elemnt, prot.nNumAtom, prot.nNumClut);
		// 力場情報データの作成を行う
//		CreateDataTop(elemnt, prot.nNumAtom);

		// タンパク質中の総原子数の計算
		prot.nNumAtom += elemnt->dataofthisresdiue.nNumAtom;
		// タンパク質中の総剛体数の計算
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
//			// 非結合性相互作用しない組
////			prot.atm[nNumAtom].not_interacting[];
//		}
//	}


}

// 座標データの作成を行う関数
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
			// 次の原子の座標を足す
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

// 剛体データの作成を行う関数
int CreateDataClust(elementoftable *element,
					int nNumAtomTotal,
					int nNumClutTotal)
{

	int nNumClust;
	int nNumBranch;
	int numChildclust;

	// 剛体原点の原子番号
	if (nNumClutTotal == 0)
	{
		prot.clt[0].origin = 1;
		prot.clt[0].nNumBranch
		= element->dataofthisresdiue.clt[0].nNumBranch;
		// 剛体が終点かどうか
		prot.clt[0].termflag
		= element->dataofthisresdiue.clt[0].termflag;
		// 剛体の原子数
		prot.clt[0].nNumAtom
		= element->dataofthisresdiue.clt[0].nNumAtom;
		// 剛体の自由度
		prot.clt[0].nDegOfFoo
		= 3;
		// 剛体終点の原子番号
		prot.clt[nNumClust+nNumClutTotal].term[nNumBranch]
		= 1;
	}

	for (nNumClust = 0;
		 nNumClust < element->dataofthisresdiue.nNumClut;
		 ++nNumClust)
	{
		if (element->dataofthisresdiue.clt[nNumClust].origin != -1)
		{
			// 剛体原点の原子番号
			prot.clt[nNumClust+nNumClutTotal].origin
				= element->dataofthisresdiue.clt[nNumClust].origin + nNumAtomTotal;
			// 剛体の枝の数
			prot.clt[nNumClust+nNumClutTotal].nNumBranch
				= element->dataofthisresdiue.clt[nNumClust].nNumBranch;
			// 剛体が終点かどうか
			prot.clt[nNumClust+nNumClutTotal].termflag
				= element->dataofthisresdiue.clt[nNumClust].termflag;
			// 剛体の原子数
			prot.clt[nNumClust+nNumClutTotal].nNumAtom
				= element->dataofthisresdiue.clt[nNumClust].nNumAtom;
			// 剛体の自由度
			prot.clt[nNumClust+nNumClutTotal].nDegOfFoo
				= 3;

			for (nNumBranch=0;
				 nNumBranch<prot.clt[nNumClust].nNumBranch;
					++nNumBranch)
			{
				// 剛体終点の原子番号
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
			// "親"は OK なので、"子"を作成
//			numChildclust = elemnt->dataofthisresdiue->clt[nNumClut+nNumClutTotal].nNumChildClut;
//			elemnt->dataofthisresdiue->clt[nNumClut+nNumClutTotal].nNumChildClut[nNumBranch] = numChildclust;
		}
	}

}

// 力場情報データの作成を行う関数
int CreateDataTop(elementoftable *element, int nNumAtomTotal)
{
	int nNumAtom;
	int nNumClust;
	int nNumDihed;
	int nNumAtomoforigin;
	int nNumClutTotal;

	// 原子種を数え上げ
	TotalNumAtomType = CountAtomID();

	// 二面角データの構築//////////////////////////////////////////////////////
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

	// 原子データの構築//////////////////////////////////////////////////////
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

// 原子種の数え上げを行う関数
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
	// n+1 番目の原子
	for (i=0;i<prot.atm[nNumAtom].num_bond_interacting;++i)
	{
		++num;
		prot.atm[nNumAtom].not_interacting[num]
		=prot.atm[nNumAtom].bond_interacting[i];
		atomnum[0]=prot.atm[nNumAtom].bond_interacting[i];
		// n+2 番目の原子
		for (i=0;i<prot.atm[atomnum[0]].num_bond_interacting;++i)
		{
			++num;
			prot.atm[nNumAtom].not_interacting[num]
			=prot.atm[nNumAtom].bond_interacting[i];
			atomnum[1]=prot.atm[nNumAtom].bond_interacting[i];
			// n+3 番目の原子
			for (i=0;i<prot.atm[atomnum[1]].num_bond_interacting;++i)
			{
				++num;
				prot.atm[nNumAtom].not_interacting[num]
				=prot.atm[nNumAtom].bond_interacting[i];
				atomnum[2]=prot.atm[nNumAtom].bond_interacting[i];
				// n+4 番目の原子
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
