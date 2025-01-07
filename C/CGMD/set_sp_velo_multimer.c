#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABA_multimer.h"
#include "gener.h"
#include "physics.h"
#include "MD.h"

// spatial velocity の設定の補助を行う関数
void sub_sub_set_sp_velo(int nNumClut,
	                     int nNumClutminusone);

// spatial velocity の設定を行う関数
void set_sp_velo(int nNumClut, int nNumClutOrigBranch){
  int i;
  int nNumClutOfParent;

  nNumClutOfParent = clust[nNumClut].nNumClutOfParent-1;

  // クラスタが末端のとき
  if (nNumClut == 0)	{
    if (TermMoveMode==ON) {
      for(i=0;i<6;++i) {
	clust[0].sp_velo[i] = vel_Term[i];
      }
    }
    else if (ABI[0].hingmat==6) {
      for(i=0;i<6;++i) {
	clust[0].sp_velo[i] = clust[0].ddihedang_six[i];
      }
    }
    else {
      for(i=0;i<6;++i) {
	clust[0].sp_velo[i] = 0.0;
      }
    }
  }
  else  {
    sub_set_sp_velo(nNumClut, nNumClutOfParent);
  }
}

// spatial velocity の設定の補助を行う関数_1
void sub_set_sp_velo(int nNumClut, int nNumClutminusone) {
  int i,j,k;
  int alpha;
  int alpha2;
  int nNumFreedom;

  FILE *test;
  // spatial velocity の初期化を行う
  for(i=0;i<6;++i) {
    clust[nNumClut].sp_velo[i] = 0.0;
  }

  // spatial velocity の設定を行う_1
  for(i=0;i<6;++i) {
    for(j=0;j<6;++j) {
      // 変換行列の転置行列を乗じる
      clust[nNumClut].sp_velo[i]
	+=  clust[nNumClut].TransMatrix[0][j][i]*clust[nNumClutminusone].sp_velo[j];
    }
  }

  // 自由度の選択
  nNumFreedom = ABI[nNumClut].hingmat-1;

  if (ABI[nNumClut].hingmat!=6)
    // spatial velocity の設定を行う_2
    clust[nNumClut].sp_velo[nNumFreedom]+= clust[nNumClut].ddihedang[0];
  else {
    for(i=0;i<6;++i) {
      clust[nNumClut].sp_velo[i]+= clust[nNumClut].ddihedang_six[i];
    }
  }

//	test = fopen("sp_velo.out", "a");
//	fprintf(test," %e \n",  clust[nNumClut].ddihedang[0]/**1.0e20*/);
//	fprintf(test,"sp_velo[%d] = ", nNumClut);
//	for (i=0;i<6;++i)
//	{
//		fprintf(test," %e \n",  clust[nNumClut].sp_velo[i]/**1.0e20*/);
//	}
//	fprintf(test, "\n");
//	fclose(test);

}



