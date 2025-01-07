
#define _GNU_SOURCE  // 2014-08-13

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "ABA.h" // 2014-06-18
#include "ABAb.h"  // 2014-06-18
#include "quaternion.h"
#include "EF.h"
#include "LA.h" // 2014-06-29
//#include "f2c.h"     // 2014-08-13
//#include "clapack.h" // 2014-08-13

#include "mymath.h" // 2014-08-13

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
             ABA_update_quaternion(deltaq[nNumClut],crd,clt,               // 2014-08-13 CHECK
			    nNumClut,num,nNumClutLast,nNumParent);  // 2014-08-13 CHECK
    }
    // In Main Chain
    else {
            ABA_update_quaternion(deltaq[nNumClut],crd,clt,                       // 2014-08-13 CHECK
      			    nNumClut,numatom-clt[nNumClut].origin_atom_a+1, // 2014-08-13 CHECK
      			    numclut-1,nNumParent);                          // 2014-08-13 CHECK
      
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
  double q_2014_08_13[4]; //q[4]; 2014-08-13
  double coord_now[4],coord_rotated[4];
  double dbmat[3][3],temp[3][3];
  double f; // 2014-08-13
  
  // 回転軸の単位v
  // 次のところは、近い将来直す。
  // 一般的には、terminal_atom_a[0]では、まずい。
  nNumAtomP = clt[nNumCltPt].terminal_atom_a[0]-1;
  nNumAtomO = clt[nNumClt].origin_atom_a-1;

  for(i=0;i<3;++i) axis_of_rotation_unit[i] = crd[nNumAtomP*3+i]-crd[nNumAtomO*3+i];
  length = 0.0;
  for(i=0;i<3;++i) length += axis_of_rotation_unit[i]*axis_of_rotation_unit[i]; 
  length = sqrt(length); // 2014-08-13 CHECK
  for(i=0;i<3;++i) axis_of_rotation_unit[i] = axis_of_rotation_unit[i]/length;

  //  q[0]=cos(delta_dihed*0.5); // 2014-08-13 CHECK
  //  printf("delta_dihed=%e %p\n",delta_dihed,&delta_dihed);  // 2014-09-05
  q_2014_08_13[0]=cos(delta_dihed*0.5); // 2014-08-13 CHECK
  //  q[0]=cos(axis_of_rotation_unit[i]); // 2014-08-13
  //  f=(double)(delta_dihed*0.5); // 2014-08-13
  //  q[0]=cos(f); // 2014-08-13 
  //  f=delta_dihed*0.5; // 2014-08-13
  //  cos(f+0.1); // 2014-08-13
  //  for(i=0;i<3;++i) q[i+1] = axis_of_rotation_unit[i]*sin(delta_dihed*0.5); // 2014-08-13 CHECK
  for(i=0;i<3;++i) q_2014_08_13[i+1] = axis_of_rotation_unit[i]*sin(delta_dihed*0.5); // 2014-08-13 CHECK

  coord_now[0]=0.0;
  for (nNumAtom=0;nNumAtom<nNumAtomALL;++nNumAtom) {
    for (i=0;i<3;++i) coord_now[i+1] = crd[(nNumAtom+nNumAtomO)*3+i]-crd[nNumAtomO*3+i];
    //    quaternion_rotation(q,coord_now,coord_rotated); // 2014-08-13
    quaternion_rotation(q_2014_08_13,coord_now,coord_rotated);
    for (i=0;i<3;++i) crd[(nNumAtom+nNumAtomO)*3+i]= coord_rotated[i+1]+crd[nNumAtomO*3+i];
  }

  // 2014-08-13
  /************************************************************/
  /* dbmat[0][0] = q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];   */
  /* dbmat[1][1] = q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];   */
  /* dbmat[2][2] = q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];   */
  /* dbmat[0][1] = 2.0*(q[1]*q[2]-q[0]*q[3]);		      */
  /* dbmat[1][0] = 2.0*(q[2]*q[1]+q[0]*q[3]);		      */
  /* dbmat[0][2] = 2.0*(q[1]*q[3]+q[0]*q[2]);		      */
  /* dbmat[2][0] = 2.0*(q[3]*q[1]-q[0]*q[2]);		      */
  /* dbmat[1][2] = 2.0*(q[2]*q[3]-q[0]*q[1]);		      */
  /* dbmat[2][1] = 2.0*(q[3]*q[2]+q[0]*q[1]);		      */
  /************************************************************/
  // 2014-08-13

  // 2014-08-13
  dbmat[0][0] = q_2014_08_13[0]*q_2014_08_13[0]+q_2014_08_13[1]*q_2014_08_13[1]-q_2014_08_13[2]*q_2014_08_13[2]-q_2014_08_13[3]*q_2014_08_13[3];
  dbmat[1][1] = q_2014_08_13[0]*q_2014_08_13[0]-q_2014_08_13[1]*q_2014_08_13[1]+q_2014_08_13[2]*q_2014_08_13[2]-q_2014_08_13[3]*q_2014_08_13[3];
  dbmat[2][2] = q_2014_08_13[0]*q_2014_08_13[0]-q_2014_08_13[1]*q_2014_08_13[1]-q_2014_08_13[2]*q_2014_08_13[2]+q_2014_08_13[3]*q_2014_08_13[3];  
  dbmat[0][1] = 2.0*(q_2014_08_13[1]*q_2014_08_13[2]-q_2014_08_13[0]*q_2014_08_13[3]);
  dbmat[1][0] = 2.0*(q_2014_08_13[2]*q_2014_08_13[1]+q_2014_08_13[0]*q_2014_08_13[3]);
  dbmat[0][2] = 2.0*(q_2014_08_13[1]*q_2014_08_13[3]+q_2014_08_13[0]*q_2014_08_13[2]);
  dbmat[2][0] = 2.0*(q_2014_08_13[3]*q_2014_08_13[1]-q_2014_08_13[0]*q_2014_08_13[2]);
  dbmat[1][2] = 2.0*(q_2014_08_13[2]*q_2014_08_13[3]-q_2014_08_13[0]*q_2014_08_13[1]);
  dbmat[2][1] = 2.0*(q_2014_08_13[3]*q_2014_08_13[2]+q_2014_08_13[0]*q_2014_08_13[1]);
  // 2014-08-13

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
void ABA_update_Term(double *crd,
		     double *delta, // 2014-08-13
		     //		     double delta[6], // 2014-08-13
		     int numatom, CLT *clt, int numclut, 
		     double *trans_A_to_CN_terminal/*[3][3]*/,double *l_Term/*[3]*/) {
  int i,j,k,l;
  int nNumAtom,nNumAtomO,nNumClt2;
  double axis_of_rotation_unit[3][3];
  double q[3][4];
  double coord_now[4],coord_rotated[4];
  double *dbmat,temp[3][3],*invmat;
  double f_temp=0.0,f; // 2014-08-13

  /******************************************************************/
  /* double *mattemp;                                 // 2014-08-13 */
  /* static long int m,n,lda,info,piv[500],lwork=500; // 2014-08-13 */
  /* static double work[500];                         // 2014-08-13 */
  /******************************************************************/

  //  dbmat=(double *)gcemalloc(sizeof(double)*3*3); // 2014-07-22
  //  dbmat=(double *)emalloc(sizeof(double)*3*3); // 2014-07-22 // 2014-09-05
  dbmat=(double *)calloc(3*3,sizeof(double)); // 2014-07-22 // 2014-09-05
  //  invmat=(double *)gcemalloc(sizeof(double)*3*3); // 2014-07-22
  //  invmat=(double *)emalloc(sizeof(double)*3*3); // 2014-07-22 // 2014-09-05
  invmat=(double *)calloc(3*3,sizeof(double)); // 2014-07-22 // 2014-09-05

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

  //  cos(0.0); // 2014-08-13

  for(i=0;i<3;++i) { // 2014-09-05
    //    q[i][0]=cos(delta[i]*0.5); // 2014-08-13
    //    cos(2.0); // 2014-08-13
    //    f_temp=axis_of_rotation_unit[2][2];  // 2014-08-13
    //    cos((double)delta[i]); // 2014-08-13
    //    printf("delta[%d]=%e %p\n",i,delta[i],&delta[i]);  // 2014-09-05
    f_temp=(double)delta[i]*0.5;       // 2014-08-13
    //    printf("f_temp=%e %p\n",f_temp,&f_temp);         // 2014-09-05
    //    printf("delta[%d]=%e %p\n",i,delta[i],&delta[i]);  // 2014-09-05

    //    f_temp=i; // 2014-08-13
    //    f_temp=(double)delta[0]*0.5;       // 2014-08-13
    q[i][0]=cos(f_temp); // 2014-08-13
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

    //    invm2(dbmat,invmat,3);                     // 2014-08-13
    //    invm(dbmat,invmat,3);                      // 2014-08-13
    //    invm3(dbmat,invmat);                             // 2014-08-13
    invm4(dbmat,invmat);                             // 2014-08-13

    /******************************************************************/
    /* m = 3;                                           // 2014-08-13 */
    /* n = 3;                                           // 2014-08-13 */
    /* lda=3;                                           // 2014-08-13 */
    /* 								      */
    /* mattemp=(double *)emalloc(sizeof(double)*m*n);   // 2014-08-13 */
    /* mtrans(dbmat,mattemp,m);                         // 2014-08-13 */
    /* dgetrf_(&m,&n,mattemp,&lda,piv,&info);           // 2014-08-13 */
    /* if (info!=0) exit(1);                            // 2014-08-13 */
    /* dgetri_(&n,mattemp,&lda,piv,work,&lwork,&info);  // 2014-08-13 */
    /* if (info!=0) exit(1);                            // 2014-08-13 */
    /* mtrans(mattemp,invmat,m);                        // 2014-08-13 */
    /* 								      */
    /* free(mattemp);                                   // 2014-08-13 */
    /******************************************************************/

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

  free(dbmat);  // 2014-07-22
  free(invmat); // 2014-07-22
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
