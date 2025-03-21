
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "ABA.h" // 2014-06-18
#include "ABAb.h"  // 2014-06-18
#include "quaternion.h"
#include "EF.h"
#include "LA.h" // 2014-06-29

void ABA_update(CLT *clt,double *crd,double *deltaq,int numclut,int numatom) {
  int i;
  int nNumClut;
  int nNumAtomLoca;
  int nNumParent;
  int nNumClut2;
  int nNumClutdummy;
  int flag,num;
  int nNumClutLast;

  for (nNumClut=0;nNumClut<numclut;++nNumClut) {
    nNumParent = clt[nNumClut].nNumClutOfParent-1;
    if (nNumClut==0) {
      ;
    }
    // In Side Chain
    else if (clt[nNumClut].join > 0) {
      nNumAtomLoca = 0;
      num=0;
      for (nNumClutdummy=nNumClut;clt[nNumClutdummy].join!=clt[nNumClut].join-1;++nNumClutdummy) {
	num+=clt[nNumClutdummy].num_atom_clust;
	nNumClutLast = nNumClutdummy;
      }
      ABA_update_quaternion(deltaq[nNumClut],crd,clt,
			    nNumClut,num,nNumClutLast,nNumParent);
    }
    // In Main Chain
    else {
      ABA_update_quaternion(deltaq[nNumClut],crd,clt,
			    nNumClut,numatom-clt[nNumClut].origin_atom_a+1,
			    numclut-1,nNumParent);
    }
  }
}

void ABA_update_quaternion(double delta_dihed,double *crd,CLT *clt,
			   int nNumClt,int nNumAtomALL, int nNumClutLast, int nNumCltPt) {
  int i,j,k,l;
  int nNumClt2,nNumClut2;
  int nNumAtom,nNumAtomO,nNumAtomP;
  double axis_of_rotation_unit[3];
  double length;
  double q[4];
  double coord_now[4],coord_rotated[4];
  double dbmat[3][3],temp[3][3];
  
  // 回転軸の単位v
  // 次のところは、近い将来直す。
  // 一般的には、terminal_atom_a[0]では、まずい。
  nNumAtomP = clt[nNumCltPt].terminal_atom_a[0]-1;
  nNumAtomO = clt[nNumClt].origin_atom_a-1;

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

// 2014-06-18
//void ABA_update_Term(double *crd,double *delta,int numatom) {
//  int i,j,k;
//
//  for (i=0;i<numatom;++i) {
//    for (j=0;j<3;++j) {
//      crd[i*3+j]+=delta[3+j];
//    }
//  }
//
//}
// 2014-06-18

// 2014-06-18
void ABA_update_Term(double *crd,double *delta, int numatom, CLT *clt, int numclut, 
		     double *trans_A_to_CN_terminal/*[3][3]*/,double *l_Term/*[3]*/) {
  int i,j,k,l;
  int nNumAtom,nNumAtomO,nNumClt2;
  double axis_of_rotation_unit[3][3];
  double q[3][4];
  double coord_now[4],coord_rotated[4];
  double *dbmat,temp[3][3],*invmat;

  dbmat=(double *)gcemalloc(sizeof(double)*3*3);
  invmat=(double *)gcemalloc(sizeof(double)*3*3);

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      ;
      //      crd[i*3+j]+=delta[3+j];
    }
  }

  axis_of_rotation_unit[0][0]=1.0/*trans_A_to_CN_terminal[0*3+0]*/;
  axis_of_rotation_unit[0][1]=0.0/*trans_A_to_CN_terminal[0*3+1]*/;
  axis_of_rotation_unit[0][2]=0.0/*trans_A_to_CN_terminal[0*3+2]*/;
			                                  	     
  axis_of_rotation_unit[1][0]=0.0/*trans_A_to_CN_terminal[1*3+0]*/;
  axis_of_rotation_unit[1][1]=1.0/*trans_A_to_CN_terminal[1*3+1]*/;
  axis_of_rotation_unit[1][2]=0.0/*trans_A_to_CN_terminal[1*3+2]*/;
			                                  	     
  axis_of_rotation_unit[2][0]=0.0/*trans_A_to_CN_terminal[2*3+0]*/;
  axis_of_rotation_unit[2][1]=0.0/*trans_A_to_CN_terminal[2*3+1]*/;
  axis_of_rotation_unit[2][2]=1.0/*trans_A_to_CN_terminal[2*3+2]*/;

  for(i=0;i<3;++i) {
    q[i][0]=cos(delta[i]*0.5);
    for(j=0;j<3;++j) q[i][j+1] = axis_of_rotation_unit[i][j]*sin(delta[i]*0.5);
  }

  //  nNumAtomO = clt[0].origin_atom_a-1;
  for (i=0;i<3;++i) l_Term[i]+=delta[i+3];

  for (nNumAtom=0;nNumAtom<numatom;++nNumAtom) {
    coord_now[0]=0.0;
    for (i=0;i<3;++i) coord_now[i+1] = crd[(nNumAtom)*3+i]/*-l_Term[i]-crd[nNumAtomO*3+i]*/;
    for (i=0;i<3;++i) quaternion_rotation(q[i],coord_now,coord_rotated);
    //    for (i=0;i<3;++i) crd[(nNumAtom)*3+i]= coord_rotated[i+1]/*+l_Term[i]+crd[nNumAtomO*3+i]*/;
  }

  for (i=0;i<3;++i) {
    dbmat[0*3+0] = q[i][0]*q[i][0]+q[i][1]*q[i][1]-q[i][2]*q[i][2]-q[i][3]*q[i][3];
    dbmat[1*3+1] = q[i][0]*q[i][0]-q[i][1]*q[i][1]+q[i][2]*q[i][2]-q[i][3]*q[i][3];
    dbmat[2*3+2] = q[i][0]*q[i][0]-q[i][1]*q[i][1]-q[i][2]*q[i][2]+q[i][3]*q[i][3];  
    dbmat[0*3+1] = 2.0*(q[i][1]*q[i][2]-q[i][0]*q[i][3]);
    dbmat[1*3+0] = 2.0*(q[i][2]*q[i][1]+q[i][0]*q[i][3]);
    dbmat[0*3+2] = 2.0*(q[i][1]*q[i][3]+q[i][0]*q[i][2]);
    dbmat[2*3+0] = 2.0*(q[i][3]*q[i][1]-q[i][0]*q[i][2]);
    dbmat[1*3+2] = 2.0*(q[i][2]*q[i][3]-q[i][0]*q[i][1]);
    dbmat[2*3+1] = 2.0*(q[i][3]*q[i][2]+q[i][0]*q[i][1]);

    for (j=0;j<3;++j) for (k=0;k<3;++k) temp[j][k] = 0.0;
    for (j=0;j<3;++j)
      for (k=0;k<3;++k)
    	for (l=0;l<3;++l)
    	  temp[j][k] += dbmat[j*3+l]*trans_A_to_CN_terminal[l*3+k/*l][k*/];

    for (j=0;j<3;++j) for (k=0;k<3;++k) trans_A_to_CN_terminal[j*3+k/*j][k*/] = temp[j][k];

    invm2(dbmat,invmat,3);

    for(nNumClt2=0;nNumClt2<numclut;++nNumClt2) {
      for (j=0;j<3;++j) for (k=0;k<3;++k) temp[j][k] = 0.0;

      for (j=0;j<3;++j) 
	for (k=0;k<3;++k) 
	  for (l=0;l<3;++l) 
	    temp[j][k] += clt[nNumClt2].trans_A_to_CN[j][l]*/*dbmat*/invmat[k*3+l];
	    //      temp[j][k] += clt[nNumClt2].trans_A_to_CN[j][l]*/*dbmat*/invmat[l*3+k];
	    //	    temp[j][k] += dbmat[j*3+l]*clt[nNumClt2].trans_A_to_CN[j][l];
	    //	    temp[j][k] += dbmat[j*3+l]*clt[nNumClt2].trans_A_to_CN[l][j];
 
      //      for (j=0;j<3;++j) for (k=0;k<3;++k) clt[nNumClt2].trans_A_to_CN[j][k] = temp[j][k];
    }
  }
}

// 2014-06-30
void ABA_update_Term2(double *crd,double *crd_Term,int numatom, double *trans_A_to_CN_terminal,double *l_Term) {
  int i,j,k,l;
  double axis_of_rotation_unit[3][3];
  double q[3][4];
  double coord_now[4],coord_rotated[4];

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) coord_rotated[j]=0.0;
    for (j=0;j<3;++j) {
      for (k=0;k<3;++k) {
	coord_rotated[j]+=trans_A_to_CN_terminal[j*3+k]*crd[i*3+k];
      }
    }
    for (j=0;j<3;++j) crd_Term[i*3+j]=coord_rotated[j];
  }

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      crd_Term[i*3+j]+=l_Term[j];
    }
  }
}
// 2014-06-30
