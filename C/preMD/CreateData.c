#include <stdio.h>
#include <stdlib.h>

#include "Vis_MD.h"

// インプットに必要なデータの作成を行う関数
// 座標データの作成を行う関数
void CreateDataCoord(elementoftable *element,int nNumAtomTotal);
// 剛体データの作成を行う関数
int CreateDataClust(elementoftable *element,int nNumAtomTotal,int nNumClutTotal);
// 力場情報データの作成を行う関数
void CreateDataTop(elementoftable *element, int nNumAtomTotal);
// 原子種の数え上げを行う関数
int CountAtomID(void);
void CkeckNonBond(int nNumAtom);
int CreateDataAtom(elementoftable *element,int nNumAtomTotal,int numTailatomonRes);
int Check(int nNumAtom, int nNumAtom2);
//int CheckDihedPair(int atoms1[4], int atoms2[4]);
void detclustdihed(int atomnum0, int numdihed);
void detttypedihed(int atomtype1, int atomtype2, int atomtype3, int atomtype4, int numdihed);
double terminatecoord[3];

// インプットに必要なデータの作成の制御を行う関数
int CreateData(void) {
	int nNumResidue;
	int AtomclustIDBEF;
	int HeadclustIDBEF=0;
	int TotalNumAtomType;
//	int nNumAtom;
//	int nNumAtom2;

	int nNumAtom;
	int nNumClut;
	int nNumDihed;

	int numTailatomonPro=0, numHeadatomonPro=0;
	int numTailatomonRes;
	int num_bond;
	int num_bond2;
	int numTailclutonPro;
	int numTailclutonRes;
	int num_branch;
//	int HAtom;

	char *nameresidue;
	elementoftable *elemnt;
	residuedata dataofthisresidue;

//	nNumAtom = 0;
//	nNumClut = 0;
//	nNumDihed = 0;

 //	for (nNumResidue=0;
 //	     nNumResidue<prot.nNumResidue;
 //	     ++nNumResidue)	{
 //	  elemnt = LURData(/*nameresidue*/nameresidue,1,dataofthisresidue);
 //	  // タンパク質中の総原子数の計算
 //	  if (nNumAtom == 0 && ADDFLAG == 1) { // to add N-H3+
 //	    nNumAtom = 2;                      //
 //	  }                                    //
 //	  nNumAtom += (elemnt->dataofthisresdiue.nNumAtom-1)/*ダミー原子の分を引いた*/;
 //	  nNumClut += elemnt->dataofthisresdiue.nNumClut;
 //	  nNumDihed += elemnt->dataofthisresdiue.nNumDihed;
 //	}

 //	prot.atm           =  (atomdata *)malloc(sizeof(atomdata)*nNumAtom);
 //	prot.atmtype       =  (atomtype *)malloc(sizeof(atomtype)*nNumAtom);
 //	prot.clt           =  (clustdata *)malloc(sizeof(clustdata)*nNumClut);
 //	prot.dihedang      =  (dihedang *)malloc(sizeof(dihedang)*nNumDihed);
 //	prot.dihedang_parm =  (dihedang_parm *)malloc(sizeof(dihedang_parm)*nNumDihed);

	prot.nNumAtom = 0;
	prot.nNumClut = 0;

	for (nNumResidue=0;
		 nNumResidue<prot.nNumResidue;
		 ++nNumResidue)
	{
		// 残基データの取得を行う
		nameresidue = prot.Sequence[nNumResidue];
		if((elemnt = LURData(nameresidue,1,dataofthisresidue))==NULL) {
		  printf("error: %s is not registered\n",nameresidue);
		  exit(1);
		}

		prot.numsequence[nNumResidue] = elemnt->dataofthisresdiue.nNumAtom;
		// 残基間の結合を行う
		// 原子
//		for (num_bond=0;;++num_bond)
//		{
//			if (prot.atm[numTailatomonPro].bond_interacting[num_bond]==0)
//			{
//				break;
//			}
//		}
//		if (nNumResidue == 0)
//		{
//			numTailatomonPro = elemnt->dataofthisresdiue.TailAtom;
//		}
		CreateDataAtom(elemnt,prot.nNumAtom,numTailatomonRes);
		if (nNumResidue > 0)
		{
			num_bond = prot.atm[numTailatomonPro-1].num_bond_interacting;
			prot.atm[numTailatomonPro-1].bond_interacting[num_bond]
			= elemnt->dataofthisresdiue.HeadAtom + prot.nNumAtom;
			prot.atm[numTailatomonPro-1].num_bond_interacting += 1;
			numHeadatomonPro = elemnt->dataofthisresdiue.HeadAtom + prot.nNumAtom;
//			numHeadatomonRes = elemnt->dataofthisresdiue.HeadAtom
			num_bond2 = elemnt->dataofthisresdiue.atm[(elemnt->dataofthisresdiue.HeadAtom)-1].num_bond_interacting;
//			num_bond = elemnt->dataofthisresdiue.HeadAtom;
//			num_bond2 = elemnt->dataofthisresdiue.atm[/*(elemnt->dataofthisresdiue.HeadAtom)*/num_bond/*-1*/].num_bond_interacting;
			prot.atm[numHeadatomonPro-1].bond_interacting[num_bond2]
			= numTailatomonPro;
			prot.atm[numHeadatomonPro-1].num_bond_interacting += 1;
		}

		// 剛体データの作成を行う
		CreateDataClust(elemnt, prot.nNumAtom, prot.nNumClut);

		// 残基間の結合を行う
		// 剛体
		if (nNumResidue > 0)
		{
			num_branch = prot.clt[numTailclutonPro].nNumBranch-1;
			prot.clt[numTailclutonPro-1].nNumChildClust[num_branch]
			= 1 + prot.nNumClut;
//			prot.clt[numTailclutonPro].nNumBranch += 1;
		}
		// 原子種を数え上げ
//		TotalNumAtomType = CountAtomID();

		// 座標データの作成を行う
		CreateDataCoord(elemnt, prot.nNumAtom);

		// 力場情報データの作成を行う
//		CreateDataTop(elemnt, prot.nNumAtom);
		numTailatomonPro = elemnt->dataofthisresdiue.TailAtom + prot.nNumAtom;
		// タンパク質中の総原子数の計算
		if (prot.nNumAtom == 0 && ADDFLAG == 1) { // to add N-H3+
		  prot.nNumAtom = 2;                      //
		}                                         //
		prot.nNumAtom += (elemnt->dataofthisresdiue.nNumAtom-1)/*ダミー原子の分を引いた*/;
		// タンパク質中の総剛体数の計算
//		if (nNumResidue < prot.nNumResidue-1)
//		{
//			prot.nNumClut += elemnt->dataofthisresdiue.nNumClut-1;
//		}
//		else
//		{
			prot.nNumClut += elemnt->dataofthisresdiue.nNumClut;
//		}
		numTailclutonPro = /*elemnt->dataofthisresdiue.nNumClut + */prot.nNumClut;
		numTailclutonRes = elemnt->dataofthisresdiue.nNumClut;
		//		num_branch = elemnt->dataofthisresdiue.atm[numTailclutonRes].nNum_Branch;

	}	 

	if (ADDFLAG == 1) {
	  prot.nNumAtom += 1;
	  prot.atm[prot.nNumAtom-1].atomID.atomnum = 0.0;
	  prot.atm[prot.nNumAtom-1].atomID.gamma = 0.0;
	  prot.atm[prot.nNumAtom-1].atomID.eata = 0.0;
	  prot.atm[prot.nNumAtom-1].atomID.e_elect = 0.0;
	  prot.atm[prot.nNumAtom-1].num_bond_interacting
	    += 1;
	  prot.atm[prot.nNumAtom-1].bond_interacting[prot.atm[prot.nNumAtom-1].num_bond_interacting]
	    = prot.nNumAtom;
	  prot.clt[prot.nNumClut-1].nNumAtom += 1;
	}
	
	prot.nNumAtomType = TotalNumAtomType;
	prot.clt[prot.nNumClut-1].nNumBranch = 0;

	// 原子種を数え上げ
	TotalNumAtomType = CountAtomID();
	// 力場情報データの作成を行う
	CreateDataTop(elemnt, TotalNumAtomType);

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

	return TotalNumAtomType;
}

// 座標データの作成を行う関数
void CreateDataCoord(elementoftable *element, int nNumAtomTotal)
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
		 nNumAtom < element->dataofthisresdiue.nNumAtom;
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
			= prot.atm[nNumAtom+nNumAtomTotal-1].coord[alpha];
	}

//	void TransLocalToAbso (int nNumAtom1, int nNumAtom2,
//                         int nSumNumAtom, double delta_dihed)

}

// 剛体データの作成を行う関数
int CreateDataAtom(elementoftable *element,
					int nNumAtomTotal,
					int numTailatomonRes)
{
	int nNumAtom;
	int nNumAtom2;
	int nNumAtom3;
	int nNumAtomMax;

	if (nNumAtomTotal==0 && ADDFLAG == 1)
	  nNumAtomMax = element->dataofthisresdiue.nNumAtom+2;
	else
	  nNumAtomMax = element->dataofthisresdiue.nNumAtom;

	for (nNumAtom=0;
	     nNumAtom<nNumAtomMax/*element->dataofthisresdiue.nNumAtom*/;
		 ++nNumAtom)
	{
		nNumAtom3 = nNumAtomTotal+nNumAtom;
		if (nNumAtom3 == 0 && ADDFLAG == 1) {
		  prot.atm[nNumAtom3].atomID.atomnum = element->dataofthisresdiue.atm[nNumAtom].atomID.atomnum;
		  prot.atm[nNumAtom3].atomID.gamma = element->dataofthisresdiue.atm[nNumAtom].atomID.gamma;
		  prot.atm[nNumAtom3].atomID.eata = element->dataofthisresdiue.atm[nNumAtom].atomID.eata;
		  prot.atm[nNumAtom3].atomID.e_elect = element->dataofthisresdiue.atm[nNumAtom].atomID.e_elect;
		  prot.atm[nNumAtom3].num_bond_interacting
		    = element->dataofthisresdiue.atm[nNumAtom].num_bond_interacting+2;
		  for (nNumAtom2=0;
		       nNumAtom2<prot.atm[/*nNumAtom*/nNumAtom3].num_bond_interacting;
		       ++nNumAtom2)  {
		    if (nNumAtom2==0 || nNumAtom2==1) {
		      prot.atm[nNumAtom3].bond_interacting[nNumAtom2] = nNumAtom2+1;
		    }
		    else {
		      prot.atm[nNumAtom3].bond_interacting[nNumAtom2]
			= element->dataofthisresdiue.atm[nNumAtom].bond_interacting[nNumAtom2]+nNumAtomTotal;
		    }
		  }
		  
		}
		else if ((nNumAtom3 == 1 || nNumAtom3 == 2) && ADDFLAG == 1) {
		  prot.atm[nNumAtom3].atomID.atomnum = 0.0;
		  prot.atm[nNumAtom3].atomID.gamma = 0.0;
		  prot.atm[nNumAtom3].atomID.eata = 0.0;
		  prot.atm[nNumAtom3].atomID.e_elect = 0.0;
		  prot.atm[nNumAtom3].num_bond_interacting
		    = 1;
		  for (nNumAtom2=0;
		       nNumAtom2<prot.atm[/*nNumAtom*/nNumAtom3].num_bond_interacting;
		       ++nNumAtom2)  {
		      prot.atm[nNumAtom3].bond_interacting[nNumAtom2]
			= element->dataofthisresdiue.clt[0].origin;
		  }
		  
		}
		else if (ADDFLAG == 1) {
		  prot.atm[nNumAtom3].atomID.atomnum = element->dataofthisresdiue.atm[nNumAtom].atomID.atomnum;
		  prot.atm[nNumAtom3].atomID.gamma = element->dataofthisresdiue.atm[nNumAtom].atomID.gamma;
		  prot.atm[nNumAtom3].atomID.eata = element->dataofthisresdiue.atm[nNumAtom].atomID.eata;
		  prot.atm[nNumAtom3].atomID.e_elect = element->dataofthisresdiue.atm[nNumAtom].atomID.e_elect;
		  prot.atm[nNumAtom3].num_bond_interacting
		    = element->dataofthisresdiue.atm[nNumAtom].num_bond_interacting;
		  for (nNumAtom2=0;
		       nNumAtom2<prot.atm[/*nNumAtom*/nNumAtom3].num_bond_interacting;
		       ++nNumAtom2)  {
		      prot.atm[nNumAtom3].bond_interacting[nNumAtom2]
			= element->dataofthisresdiue.atm[nNumAtom].bond_interacting[nNumAtom2]+nNumAtomTotal+2;
		  }

		}
		else {
		  prot.atm[nNumAtom3].atomID.atomnum = element->dataofthisresdiue.atm[nNumAtom].atomID.atomnum;
		  prot.atm[nNumAtom3].atomID.gamma = element->dataofthisresdiue.atm[nNumAtom].atomID.gamma;
		  prot.atm[nNumAtom3].atomID.eata = element->dataofthisresdiue.atm[nNumAtom].atomID.eata;
		  prot.atm[nNumAtom3].atomID.e_elect = element->dataofthisresdiue.atm[nNumAtom].atomID.e_elect;
		  prot.atm[nNumAtom3].num_bond_interacting
		    = element->dataofthisresdiue.atm[nNumAtom].num_bond_interacting;
		  for (nNumAtom2=0;
		       nNumAtom2<prot.atm[/*nNumAtom*/nNumAtom3].num_bond_interacting;
		       ++nNumAtom2)  {
		      prot.atm[nNumAtom3].bond_interacting[nNumAtom2]
			= element->dataofthisresdiue.atm[nNumAtom].bond_interacting[nNumAtom2]+nNumAtomTotal;
		  }
		}
	}
}

// 剛体データの作成を行う関数
int CreateDataClust(elementoftable *element,
					int nNumAtomTotal,
					int nNumClutTotal)
{

	int nNumClust;
	int nNumBranch;
	int numChildclust;

//	// 剛体原点の原子番号
//	if (nNumClutTotal == 0)
//	{
//		prot.clt[0].origin = 1;
//		prot.clt[0].nNumBranch
//		= element->dataofthisresdiue.clt[0].nNumBranch;
//		// 剛体が終点かどうか
//		prot.clt[0].termflag
//		= element->dataofthisresdiue.clt[0].termflag;
//		// 剛体の原子数
//		prot.clt[0].nNumAtom
//		= element->dataofthisresdiue.clt[0].nNumAtom;
//		// 剛体の自由度
//		prot.clt[0].nDegOfFoo
//		= 3;
//		// 剛体終点の原子番号
//		prot.clt[nNumClust+nNumClutTotal].term[nNumBranch]
//		= 1;
//	}
		for (nNumClust = 0;
			 nNumClust < element->dataofthisresdiue.nNumClut;
			 ++nNumClust)
		{
//		if (element->dataofthisresdiue.clt[nNumClust].origin != -1)
//		{
			// 剛体原点の原子番号
			prot.clt[nNumClust+nNumClutTotal].origin
				= element->dataofthisresdiue.clt[nNumClust].origin + nNumAtomTotal;
			if (nNumAtomTotal == 0 && nNumClust != 0 && ADDFLAG == 1) {
			  prot.clt[nNumClust+nNumClutTotal].origin += 2;
			}
			// 剛体の枝の数
			prot.clt[nNumClust+nNumClutTotal].nNumBranch
				= element->dataofthisresdiue.clt[nNumClust].nNumBranch;
			if (nNumClust == 0 && nNumClutTotal > 0)
			{
//				// 剛体の枝の数
//				prot.clt[nNumClutTotal-1].nNumBranch
//					+= 1;
			}
			// 剛体が終点かどうか
			prot.clt[nNumClust+nNumClutTotal].termflag
				= element->dataofthisresdiue.clt[nNumClust].termflag;
			if (nNumClust == 0 && nNumClutTotal > 0)
			{
			// 剛体が終点かどうか
				prot.clt[nNumClutTotal-1].termflag
				= 1;
			}
			// 剛体の原子数
			prot.clt[nNumClust+nNumClutTotal].nNumAtom
				= element->dataofthisresdiue.clt[nNumClust].nNumAtom;
			if (nNumAtomTotal == 0 && nNumClust == 0 && ADDFLAG == 1) {
			  prot.clt[nNumClust+nNumClutTotal].nNumAtom += 2;
			}
			// 剛体の自由度
			prot.clt[nNumClust+nNumClutTotal].nDegOfFoo
				= 3;

			for (nNumBranch=0;
				 nNumBranch<prot.clt[/*nNumClust*/nNumClust+nNumClutTotal].nNumBranch;
					++nNumBranch)
			{
				// 剛体終点の原子番号
				prot.clt[nNumClust+nNumClutTotal].term[nNumBranch]
				=   element->dataofthisresdiue.clt[nNumClust].term[nNumBranch]
				   + nNumAtomTotal;
				if (nNumAtomTotal == 0 && nNumClust != 0 && ADDFLAG == 1) {
				  prot.clt[nNumClust+nNumClutTotal].term[nNumBranch] += 2;
				}
			}
//		}
	}
	for (nNumClust = 0;
		 nNumClust < element->dataofthisresdiue.nNumClut;
		 ++nNumClust)
	{
		if (nNumClust == 0)
		{
//			prot.clt[nNumClutTotal-1].nNumParentClust += 1;
		}
		prot.clt[nNumClust+nNumClutTotal].nNumParentClust
		= element->dataofthisresdiue.clt[nNumClust].nNumParentClust + nNumClutTotal;
		for (nNumBranch=0;
			 nNumBranch<prot.clt[/*nNumClust*/nNumClust+nNumClutTotal].nNumBranch;
				++nNumBranch)
		{
			// "親"は OK なので、"子"を作成
			prot.clt[nNumClust+nNumClutTotal].nNumChildClust[nNumBranch]
			= element->dataofthisresdiue.clt[nNumClust].nNumChildClust[nNumBranch] + nNumClutTotal;
		}
	}

}

// 力場情報データの作成を行う関数
void CreateDataTop(elementoftable *element, int TotalNumAtomType)
{
	int nNumAtom;
	int nNumClust;
	int nNumDihed;
	int nNumAtomoforigin;
	int nNumClutTotal;

//	prot.nNumDihedType = 0;
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
		 nNumAtom<prot.nNumAtom;
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
	int MaxNumAtomType=40;
	int ThisAtomType;
	int numAtom;
	int flag = 1;
	int ON = 1;
	int OFF = 0;
	int list[20/*20*//*MAXNAT*/];
	int i,j;

	for (i=0;i<20/*20*/;++i)
	{
		list[i]=0;
	}

//	for (numAtom=0;numAtom<prot.nNumAtom;++numAtom)
//	{
//		ThisAtomType = prot.atm[numAtom].atomID.atomnum;
//		for (j=0;j<=/*Total*/MaxNumAtomType;++j)
//		{
//			if (ThisAtomType == list[j])
//			{
//				flag = OFF;
//				break;
//			}
//		}
//		if (flag == ON)
//		{
//			list[TotalNumAtomType] = ThisAtomType;
//			prot.parmA[TotalNumAtomType] = prot.atm[numAtom].atomID.gamma;
//			prot.parmB[TotalNumAtomType] = prot.atm[numAtom].atomID.eata;
//			++TotalNumAtomType;
//		}
//		flag = ON;
//	}

	for (numAtom=0;numAtom<20;++numAtom)
	{
		prot.parmA[numAtom] = 0.0;
		prot.parmB[numAtom] = 0.0;
	}

	for (numAtom=0;numAtom<prot.nNumAtom;++numAtom)
	{
		ThisAtomType = prot.atm[numAtom].atomID.atomnum;
		prot.parmA[ThisAtomType-1] = prot.atm[numAtom].atomID.gamma;
		prot.parmB[ThisAtomType-1] = prot.atm[numAtom].atomID.eata;
	}

	return /*TotalNumAtomType*/20;
}

void CkeckNonBond(int nNumAtom)
{
	int i,j,k,l,m;
	int num=0,num2=0;
	int atomnum[3];
	int atomnum2[3];
	int flag;
	int ON=0;
	int OFF=1;
	int atomtype1;
	int atomtype2;
	int CheckDihed;
	int numdihedsofar;

	prot.atm[nNumAtom].not_interacting[/*0*/num]=nNumAtom+1;
	// n+1 番目の原子
	for (i=0;i<prot.atm[nNumAtom].num_bond_interacting;++i)
	{
		flag = ON;
		if (Check(nNumAtom, prot.atm[nNumAtom].bond_interacting[i])==0)
		{
			++num;
			prot.atm[nNumAtom].not_interacting[num]
			=prot.atm[nNumAtom].bond_interacting[i];
		}
//		else
//			break;
		atomnum[0]=prot.atm[nNumAtom].bond_interacting[i]-1;
		atomnum2[0]=prot.atm[nNumAtom].bond_interacting[i]-1;
		// n+2 番目の原子
		for (j=0;j<prot.atm[atomnum[0]].num_bond_interacting;++j)
		{
			if (Check(nNumAtom, prot.atm[atomnum[0]].bond_interacting[j])==0)
			{
				++num;
				prot.atm[nNumAtom].not_interacting[num]
				=prot.atm[atomnum[0]].bond_interacting[j];
			}
//			else
//				break;
			if (atomnum[0] != prot.atm[atomnum[0]].bond_interacting[j]-1
			&&  nNumAtom   != prot.atm[atomnum[0]].bond_interacting[j]-1)
			{
				atomnum[1]=prot.atm[atomnum[0]].bond_interacting[j]-1;
				if (prot.atm[atomnum2[0]].bond_interacting[j]-1 > atomnum2[0])
				{
					atomnum2[1]=prot.atm[atomnum2[0]].bond_interacting[j]-1;
				}
				else
				{
					flag = OFF;
				}
				// n+3 番目の原子
				for (k=0;k<prot.atm[atomnum[1]].num_bond_interacting;++k)
				{
					if (Check(nNumAtom, prot.atm[atomnum[1]].bond_interacting[k])==0)
					{
						++num;
						prot.atm[nNumAtom].not_interacting[num]
						=prot.atm[atomnum[1]].bond_interacting[k];
					}
//				else
//					break;
					if (atomnum[0] != prot.atm[atomnum[1]].bond_interacting[k]-1
					 &&  nNumAtom  != prot.atm[atomnum[1]].bond_interacting[k]-1
					 && atomnum[1] != prot.atm[atomnum[1]].bond_interacting[k]-1)
					{
						atomnum[2]=prot.atm[atomnum[1]].bond_interacting[k]-1;
						if (prot.atm[atomnum2[1]].bond_interacting[k]-1 > atomnum2[1])
						{
							atomnum2[2]=prot.atm[atomnum2[1]].bond_interacting[k]-1;
						}
						else
						{
							flag = OFF;
						}
						// n+4 番目の原子
						for (l=0;l<prot.atm[atomnum[2]].num_bond_interacting;++l)
						{
							if (Check(nNumAtom, prot.atm[atomnum[2]].bond_interacting[l])==0)
							{
//								++num;
//								prot.atm[nNumAtom].not_interacting[num]
//								=prot.atm[atomnum[2]].bond_interacting[l];
							}
//					else
//						break;
							if (Check2(nNumAtom, atomnum[2]+1)==0 /*&& flag == ON*/)
							{
								prot.atm[nNumAtom].o_f_interacting[num2]
//								=prot.atm[atomnum2[2]].bond_interacting[l];
								=atomnum[2]+1;
								CheckDihed = 0;
/*								for(numdihedsofar=0;numdihedsofar<num2+prot.nNumDihedAng;++numdihedsofar)
								{
									if(CheckDihedPair(prot.dihedang[numdihedsofar].atomnum,atomnum) == 1)
									{
										CheckDihed = 1;
										break;
									}
								}
								if(CheckDihed != 1)
								{*/
									prot.dihedang[num2+prot.nNumDihedAng].atomnum[0] = nNumAtom;
									for (m=1;m<4;++m)
									{
										prot.dihedang[num2+prot.nNumDihedAng].atomnum[m] = atomnum[m-1];
									}
									detclustdihed(prot.dihedang[num2+prot.nNumDihedAng].atomnum[2],
												  num2+prot.nNumDihedAng);
									detttypedihed(prot.atm[nNumAtom].atomID.atomnum,
												  prot.atm[atomnum[0]].atomID.atomnum,
												  prot.atm[atomnum[1]].atomID.atomnum,
												  prot.atm[atomnum[2]].atomID.atomnum, 
												  num2+prot.nNumDihedAng);
									++num2;
/*								}*/
							}
						}
					}
				}
			}
		}
	}

	prot.nNumDihedAng += num2;
	prot.atm[nNumAtom].num_not_interacting=num;
	prot.atm[nNumAtom].num_o_f_interacting=num2;
}

int Check(int nNumAtom, int nNumAtom2)
{
	int i;

	for (i=0;;++i)
	{
		if (nNumAtom2 == prot.atm[nNumAtom].not_interacting[i])
		{
			return 1;
//			break;
		}
		else if (prot.atm[nNumAtom].not_interacting[i] == 0 /*&& nNumAtom > 0*/)
		{
			return 0;
			break;
		}
	}

//	return 0;
}

int Check2(int nNumAtom, int nNumAtom2)
{
	int i;

	for (i=0;;++i)
	{
		if (nNumAtom2 == prot.atm[nNumAtom].o_f_interacting[i])
		{
			return 1;
			break;
		}
		else if (prot.atm[nNumAtom].o_f_interacting[i] == 0)
		{
			break;
		}
	}

	return 0;
}

void detttypedihed(int atomtype1, int atomtype2, int atomtype3, int atomtype4, int numdihed)
{
	int DihedType;
	int numDihedType;

	for (numDihedType=0;numDihedType<20;++numDihedType)
	{
		if ((atomtype1==dataofdihed[numDihedType].atomnum[0]
		  && atomtype2==dataofdihed[numDihedType].atomnum[1]
		  && atomtype3==dataofdihed[numDihedType].atomnum[2]
		  && atomtype4==dataofdihed[numDihedType].atomnum[3])||
		    (atomtype1==dataofdihed[numDihedType].atomnum[3]
		  && atomtype2==dataofdihed[numDihedType].atomnum[2]
		  && atomtype3==dataofdihed[numDihedType].atomnum[1]
		  && atomtype4==dataofdihed[numDihedType].atomnum[0]))
		{
			prot.dihedang[numdihed].atomnum[4] = numDihedType;
			prot.dihedang_parm[numdihed].num = dataofdihed[numDihedType].dihedangp.num;
			prot.dihedang_parm[numdihed].V_2 = dataofdihed[numDihedType].dihedangp.V_2;
			prot.dihedang_parm[numdihed].theta = dataofdihed[numDihedType].dihedangp.theta;
//			prot.nNumDihedType += 1;
		}
	}
}

void detclustdihed(int atomnum0, int numdihed)
{
	int numclust;

	for (numclust=0;numclust<prot.nNumClut;++numclust)
	{
		if(atomnum0 == prot.clt[numclust].origin-1)
		{
			prot.dihedang[numdihed].atomnum[5] = numclust+1;
			break;
		}
	}

}

/*int CheckDihedPair(int atoms1[4], int atoms2[4])
{
	int i;

	for (i=0;atoms2[3-i]==atoms1[i];++i)
	{
		if(i==3)
		{
			return 1;
		}
	}

	if (i<3)
	return 0;
}*/

