
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABAb.h"
#include "quaternion.h"
#include "EF.h"

void ABAb_update(CLTb *clt,double *crd,double *deltaq,int numclut,int numatom) {
  int i;
  int nNumClut;
  int nNumAtomLoca;
  int nNumParent;
  int nNumClut2;
  int nNumClutdummy;
  int flag,num;
  int nNumClutLast;
  int nNumAtomP,nNumAtomO;
  int *num_index_terminal_atom_a;

  num_index_terminal_atom_a=(int)gcemalloc(sizeof(int)*numclut);
  for (nNumClut=0;nNumClut<numclut;++nNumClut) {
    num_index_terminal_atom_a[nNumClut]=0;
  }

  for (nNumClut=0;nNumClut<numclut;++nNumClut) {
    nNumParent = clt[nNumClut].nNumClutOfParent-1;
    if (nNumClut!=0) {
      i=num_index_terminal_atom_a[nNumParent];
      ++num_index_terminal_atom_a[nNumParent];
      nNumAtomP = clt[nNumParent].terminal_atom_a[i]-1;
      nNumAtomO = clt[nNumClut].origin_atom_a-1;
    }
    if (nNumClut==0) {
      ;
    }
    // In Side Chain
    else if (clt[nNumClut].join > 0) {
      nNumAtomLoca = 0;
      num=0;
      for (nNumClutdummy=nNumClut;
	   /*clt[nNumClutdummy].join!=clt[nNumClut].join-1*/;
	   ++nNumClutdummy) {
	num+=clt[nNumClutdummy].num_atom_clust;
	nNumClutLast = nNumClutdummy;
	if ((clt[nNumClutdummy].terminal==0) && (clt[nNumClutdummy].join==clt[nNumClut].join))
	  break;
      }
      ABAb_update_quaternion(deltaq[nNumClut],crd,clt,
			     nNumClut,num,nNumClutLast,nNumParent,
			     nNumAtomP,nNumAtomO);
    }
    // In Main Chain
    else {
      ABAb_update_quaternion(deltaq[nNumClut],crd,clt,
			     nNumClut,numatom-clt[nNumClut].origin_atom_a+1,
			     numclut-1,nNumParent,
			     nNumAtomP,nNumAtomO);
    }
  }
}

void ABAb_update_quaternion(double delta_dihed,double *crd,CLTb *clt,
			    int nNumClt,int nNumAtomALL, int nNumClutLast, int nNumCltPt,
			    int nNumAtomP, int nNumAtomO) {
  int i,j,k,l;			
  int nNumClt2,nNumClut2;
  int nNumAtom/*,nNumAtomO,nNumAtomP*/;
  double axis_of_rotation_unit[3];
  double length;
  double q[4];
  double coord_now[4],coord_rotated[4];
  double dbmat[3][3],temp[3][3];
  
  // 回転軸の単位v
  // 次のところは、近い将来直す。
  // 一般的には、terminal_atom_a[0]では、まずい。
  //  nNumAtomP = clt[nNumCltPt].terminal_atom_a[0]-1;
  //  nNumAtomO = clt[nNumClt].origin_atom_a-1;

  for(i=0;i<3;++i) axis_of_rotation_unit[i] = crd[nNumAtomP*3+i]-crd[nNumAtomO*3+i];
  length = 0.0;
  for(i=0;i<3;++i) length += axis_of_rotation_unit[i]*axis_of_rotation_unit[i]; 
  length = sqrt(length);
  for(i=0;i<3;++i) axis_of_rotation_unit[i] = axis_of_rotation_unit[i]/length;

  q[0]=cos(delta_dihed*0.5);
  for(i=0;i<3;++i) q[i+1] = axis_of_rotation_unit[i]*sin(delta_dihed*0.5);

  coord_now[0]=0.0;
  for (nNumAtom=0;nNumAtom<nNumAtomALL;++nNumAtom) {
    for (i=0;i<3;++i) coord_now[i+1] = crd[(nNumAtom+nNumAtomO)*3+i]-crd[nNumAtomO*3+i];
    quaternion_rotation(q,coord_now,coord_rotated);
    for (i=0;i<3;++i) crd[(nNumAtom+nNumAtomO)*3+i]= coord_rotated[i+1]+crd[nNumAtomO*3+i];
  }

  dbmat[0][0] = q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
  dbmat[1][1] = q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
  dbmat[2][2] = q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];  
  dbmat[0][1] = 2.0*(q[1]*q[2]-q[0]*q[3]);
  dbmat[1][0] = 2.0*(q[2]*q[1]+q[0]*q[3]);
  dbmat[0][2] = 2.0*(q[1]*q[3]+q[0]*q[2]);
  dbmat[2][0] = 2.0*(q[3]*q[1]-q[0]*q[2]);
  dbmat[1][2] = 2.0*(q[2]*q[3]-q[0]*q[1]);
  dbmat[2][1] = 2.0*(q[3]*q[2]+q[0]*q[1]);

  for(nNumClt2=nNumClt;nNumClt2<=nNumClutLast;++nNumClt2) {
    for (i=0;i<3;++i) 
      for (j=0;j<3;++j) temp[i][j] = 0.0;
    for (i=0;i<3;++i) 
      for (j=0;j<3;++j) 
	for (k=0;k<3;++k) 
	  temp[i][j] += clt[nNumClt2].trans_A_to_CN[i][k]*dbmat[j][k];
 
    for (i=0;i<3;++i) 
      for (j=0;j<3;++j) 
    	clt[nNumClt2].trans_A_to_CN[i][j] = temp[i][j];
  }
}

void ABAb_update_Term(double *crd,double *delta,int numatom) {
  int i,j,k;

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      crd[i*3+j]+=delta[3+j];
    }
  }

}

