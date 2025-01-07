#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prot.h"
#include "quaternion.h"

void sub_trans_CN_to_A_quaternion(int nNumClt, int nNumAtomALL, int nNumClutLast, int nNumCltPt);

// 局所座標系→実験室系の変換を行う関数
void trans_CN_to_A(int nNumClut, int nNumClutOrigBranch) {
  int i;
  int nNumAtomLoca;
  int nNumParent;
  int nNumClut2;
  int num_atom;
  int nNumClutdummy;
  int flag,num;
  int nNumClutLast;

  nNumParent = clust[nNumClut].nNumClutOfParent-1;

  // 0番目のクラスタではなにもしない
  if (nNumClut == 0) {
    ;
  }
  // In Side Chain
  else if (clust[nNumClut].join > 0) {
    nNumAtomLoca = 0;
    num=0;
    for (nNumClutdummy=nNumClut;clust[nNumClutdummy].join!=clust[nNumClut].join-1;++nNumClutdummy) {
      num+=clust[nNumClutdummy].num_atom_clust;
      nNumClutLast = nNumClutdummy;
    }
    // 局所座標系→実験室系の変換
       sub_trans_CN_to_A_quaternion(nNumClut,num,nNumClutLast,/*nNumAtomLoca*/nNumParent);
  }
  // In Main Chain
  else {
    // 局所座標系→実験室系の変換
    sub_trans_CN_to_A_quaternion(nNumClut,prot.num_atom-clust[nNumClut].origin_atom_a+1,prot.DOF-1,/*0*/nNumParent);
  }
	
}

// 局所座標系→実験室系の変換を行う関数
void sub_trans_CN_to_A_quaternion(int nNumClt,int nNumAtomALL, int nNumClutLast, int nNumCltPt) {
  int i,j,k,l;
  int nNumClt2,nNumClut2;
  int nNumAtom;
  int nNumAtomO;
  int nNumAtomP;
  double delta_dihed;
  double axis_of_rotation_unit[3];
  double length;
  double q[4];
  double coord_now[4],coord_rotated[4];
  double dbmat[3][3],dbmattemp[3][3],dbmattemp2[3][3],dbmat2[3][3],temp[3][3];
  
  // 現ステップでの二面角の変位
  delta_dihed = clust[nNumClt].now_deltadihedang[0];
  // 回転軸の単位v
  // 次のところは、近い将来直す。
  // 一般的には、terminal_atom_a[0]では、まずい。
  nNumAtomP = clust[nNumCltPt].terminal_atom_a[0]-1;
  nNumAtomO = clust[nNumClt].origin_atom_a-1;


  for(i=0;i<3;++i) {
    axis_of_rotation_unit[i] = prot.coord[nNumAtomP][i]-prot.coord[nNumAtomO][i];
  }
  length = 0.0;
  for(i=0;i<3;++i) {
    length += axis_of_rotation_unit[i]*axis_of_rotation_unit[i]; 
  }
  length = sqrt(length);
  for(i=0;i<3;++i) {
    axis_of_rotation_unit[i] = axis_of_rotation_unit[i]/length;
  }

  q[0]=cos(delta_dihed*0.5);
  for(i=0;i<3;++i) {
    q[i+1] = axis_of_rotation_unit[i]*sin(delta_dihed*0.5);
  }

  coord_now[0]=0.0;
  for (nNumAtom=0;nNumAtom<nNumAtomALL;++nNumAtom) {
    for (i=0;i<3;++i) {
      coord_now[i+1] =  
	            prot.coord[nNumAtom+nNumAtomO][i]
                   -prot.coord[nNumAtomO][i]
	;
    }

    quaternion_rotation(q,coord_now,coord_rotated);
    for (i=0;i<3;++i) {
      prot.coord[nNumAtom+nNumAtomO][i]= coord_rotated[i+1]  
                                       +prot.coord[nNumAtomO][i]
	;
    }
    

  }

  // 次のステップの座標での二面角の設定
  clust[nNumClt].dihedang[0] += delta_dihed;
  
    
  dbmat[0][0] = q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
  dbmat[1][1] = q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
  dbmat[2][2] = q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
  
  dbmat[0][1] = 2.0*(q[1]*q[2]-q[0]*q[3]);
  dbmat[1][0] = 2.0*(q[2]*q[1]+q[0]*q[3]);
  
  dbmat[0][2] = 2.0*(q[1]*q[3]+q[0]*q[2]);
  dbmat[2][0] = 2.0*(q[3]*q[1]-q[0]*q[2]);
  
  dbmat[1][2] = 2.0*(q[2]*q[3]-q[0]*q[1]);
  dbmat[2][1] = 2.0*(q[3]*q[2]+q[0]*q[1]);
    
  // 次のステップの座標での局所座標系の計算
  for(nNumClt2=nNumClt;nNumClt2<=nNumClutLast; ++nNumClt2) {
    for (i=0;i<3;++i) {
      for (j=0;j<3;++j) {
	temp[i][j] = 0.0;
      }
    }
    for (i=0;i<3;++i) {
      for (j=0;j<3;++j) {
	for (k=0;k<3;++k) {
	  temp[i][j] += clust[nNumClt2].trans_A_to_CN[0][i][k]*dbmat[j][k];
	}
      }
    }
    for (i=0;i<3;++i) {
      for (j=0;j<3;++j) {
    	clust[nNumClt2].trans_A_to_CN[0][i][j] = temp[i][j];
      }
    }
  }
}
