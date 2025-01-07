#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "Vis_MD.h"

// 関数
void atomscanf( FILE *dataf, /* residuedata dataofthisresidue */int NumRes);
void clustscanf(FILE *dataf, /* residuedata dataofthisresidue */int NumRes);
void dihedscanf(FILE *dataf, elementoftable *data, residuedata dataofthisresidue);

// データの取得を行う関数
int PickData(void) {
  int nNumRes;
  char	nameofthisresidue[MAXNRESIDUE][4];
  int	flag;
  FILE *datafile;
  
  // 標準残基データベースを開く
  datafile = efopen("stdres.dat","r");
  
  // 標準残基データベースの読み込み
  for (nNumRes=0 ; nNumRes < NRESIDUE ; ++nNumRes){
    // 残基名を読む
    fscanf(datafile, "%s", &(nameofthisresidue[nNumRes]));
    nameofthisresidue[nNumRes][3] = '\0';
    // 残基中の原子数を読む
    fscanf(datafile, "%d", &(dataofresidueontable[nNumRes].nNumAtom));
    // 残基中の原子のデータを読む
    atomscanf(datafile,/* dataofresidueontable[nNumRes] */nNumRes);
    // 残基中の剛体数を読む
    fscanf(datafile, "%d", &(dataofresidueontable[nNumRes].nNumClut));
    // 残基中の剛体のデータを読む
    clustscanf(datafile,/* dataofresidueontable[nNumRes]*/nNumRes);
    // テーブル中にデータを登録
    LURData(nameofthisresidue[nNumRes],0,dataofresidueontable[nNumRes]);
    fscanf(datafile, "%d", &flag);

    if (++NumResInStdResData > NRESIDUE) {
      printf("too many residues in stdres.dat !!!\n");
      break;
    }
    
    if (flag == 0)
      break;
  }

  // 標準残基データベースを閉じる
  fclose(datafile);
}

// 標準残基中の原子データの読み込みを行う関数
void atomscanf(FILE *dataf, /* residuedata dataofthisresidue */int nNumRes ) {
  int alpha;
  int nNumAtom;
  int nNumAtom2;

  int nNum_not_interacting;
  int nNum_not_interacting_kminusone;
  int nNum_not_interacting_kplusone;
  
  int ThisnNumAtom;
  
  int nNumAtomOfResidue;
  
  double dData;
  double e_elect;
  double eata;
  double gamma;
  int nData;
  char sData;

  nNumAtomOfResidue = dataofresidueontable[nNumRes].nNumAtom;

  // 通常原子のデータの読み込み
  for (nNumAtom = 0; nNumAtom < nNumAtomOfResidue; ++nNumAtom) {
    // 原子番号の読み込み
    fscanf(dataf, "%d", &ThisnNumAtom);
    // ダミー原子の場合、原子情報はナシ
    if (ThisnNumAtom != nNumAtomOfResidue) 
      // 原子名の読み込み
      fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].atm[nNumAtom].atomID.atomnum));
    
    // head or tail 原子の読み込み
    fscanf(dataf, "%d", &nData);
    if (nData == 1)
      dataofresidueontable[nNumRes].HeadAtom = nNumAtom+1;
    else if (nData == -1)
	dataofresidueontable[nNumRes].TailAtom = nNumAtom+1;

    // 座標データの読み込み
    for (alpha = 0; alpha < 3; ++alpha)
      fscanf(dataf, "%lf", &(dataofresidueontable[nNumRes].atm[nNumAtom].coord[alpha]));

    // ダミー原子の場合、原子情報はナシ
    if (ThisnNumAtom != nNumAtomOfResidue) {
      // 電荷の読み込み
      fscanf(dataf, "%lf", &(dataofresidueontable[nNumRes].atm[nNumAtom].atomID.e_elect));
      // ガンマの読み込み
      fscanf(dataf, "%lf", &(dataofresidueontable[nNumRes].atm[nNumAtom].atomID.gamma));
      // イータの読み込み
      fscanf(dataf, "%lf", &(dataofresidueontable[nNumRes].atm[nNumAtom].atomID.eata));      
      fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].atm[nNumAtom].num_bond_interacting));
      // 結合の原子の読み込み
      for (nNumAtom2=0;nNumAtom2<dataofresidueontable[nNumRes].atm[nNumAtom].num_bond_interacting;++nNumAtom2) {
	fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].atm[nNumAtom].bond_interacting[nNumAtom2]));
      }
    }
  }
}


// 標準残基中の剛体データの読み込みを行う関数
void clustscanf(FILE *dataf, /* residuedata dataofthisresidue */int nNumRes ) {
  int nNumClt;
  int nNumBh;
  int nNumCltOfResidue;
  
  nNumCltOfResidue = dataofresidueontable[nNumRes].nNumClut;
  // 剛体データの読み込み
  for (nNumClt = 0; nNumClt < nNumCltOfResidue; ++nNumClt) {
    // 剛体中の原子数の読み込み
    fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].clt[nNumClt].nNumAtom));    
    // 剛体中の枝の数の読み込み
    fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].clt[nNumClt].nNumBranch));        
    // 剛体の末端情報の読み込み
    fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].clt[nNumClt].termflag));    
    // 剛体中の原点情報の読み込み
    fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].clt[nNumClt].origin));
    // 剛体中の終点情報の読み込み
    for (nNumBh = 0; nNumBh < dataofresidueontable[nNumRes].clt[nNumClt].nNumBranch; ++nNumBh)
      fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].clt[nNumClt].term[nNumBh]));	
    // "親"の剛体の番号の読み込み
    fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].clt[nNumClt].nNumParentClust));    
    // "子"の剛体の番号の読み込み
    for (nNumBh = 0; nNumBh < dataofresidueontable[nNumRes].clt[nNumClt].nNumBranch; ++nNumBh)
      fscanf(dataf, "%d", &(dataofresidueontable[nNumRes].clt[nNumClt].nNumChildClust[nNumBh]));
  }
}

// 標準残基中の二面角データの読み込みを行う関数
void dihedscanf(FILE *dataf,elementoftable *data, residuedata dataofthisresidue) {
  int nNumDihed;
  int nNumAtm;
  double dData;
  int nData;
  int nNumDihedOfResidue;

  fscanf(dataf, "%d", &(prot.nNumDihedType));

  // 二面角データの読み込み
  for (nNumDihed = 0; nNumDihed < nNumDihedOfResidue; ++nNumDihed) {
    fscanf(dataf, "%d", &nData);    
    // 二面角相互作用のV_2の読み込み
    fscanf(dataf, "%lf", &(dataofdihed[nNumDihed].dihedangp.V_2));
    // 二面角相互作用の多重度の読み込み
    fscanf(dataf, "%d", &(dataofdihed[nNumDihed].dihedangp.num));
    // 二面角相互作用の位相の読み込み
    fscanf(dataf, "%lf", &(dataofdihed[nNumDihed].dihedangp.theta));
    for (nNumAtm=0;nNumAtm<4;++nNumAtm)	{
      fscanf(dataf, "%d", &(dataofdihed[nNumDihed].atomnum[nNumAtm]));
    }
  }
}

