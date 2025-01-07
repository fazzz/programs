#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gener.h"
#include "MD.h"
#include "ABA_multimer.h"
#include "quaternion.h"

#include "EF.h"

//#include "QUA.h"

//void mvult(double m1[6][6], double v[6], double mv[6], int n);
//void mmult(double m1[6][6], double m2[6][6], double m1m2[6][6], int n);

//void sub_trans_CN_to_A_term(void);
//void sub_trans_CN_to_A_term(double q[4],double omg[3]);
void sub_trans_CN_to_A_term(double q_Term[4]);
void sub_trans_CN_to_A_quaternion(int nNumClt, int nNumAtomALL, int nNumClutLast, int nNumCltPt);

// 局所座標系→実験室系の変換を行う関数
void trans_CN_to_A(int nNumClut, int nNumClutOrigBranch, double q_Term[4]) {
  int i;
  int nNumAtomLoca;
  int nNumParent;
  int nNumClut2;
  int num_atom;
  int nNumClutdummy;
  int flag,num;
  int nNumClutLast;

  int nNumClutALL,nNumAtomALL;

  nNumParent = clust[nNumClut].nNumClutOfParent-1;


  // 0番目のクラスタではなにもしない
  if (nNumClut == 0) {
    if (TermMoveMode == ON ) {
      sub_trans_CN_to_A_term(q_Term);
    }
    else if ( ABI[0].hingmat==6 ) {
      nNumAtomALL=clust[nNumClut].num_atom_clust;
      nNumClutALL=0;
      sub_trans_CN_to_A_quaternion_six(nNumClut,nNumAtomALL,nNumClutALL,/*0*/nNumParent);
    }
    else {
      ;
    }

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
    if (ABI[nNumClut].hingmat!=6)
      sub_trans_CN_to_A_quaternion(nNumClut,num,nNumClutLast,/*nNumAtomLoca*/nNumParent);
    else {
      //      sub_trans_CN_to_A_quaternion_six(nNumClut,num,nNumClutLast,/*nNumAtomLoca*/nNumParent);
      ;
       //    sub_trans_CN_to_A(nNumClut,num,nNumAtomLoca);
    } 
  }
  // In Main Chain
  else {
    // 局所座標系→実験室系の変換
    nNumAtomALL=clust[nNumClut].num_atom_clust;
    for (nNumClutdummy=0;
	 ABI[nNumClut+1+nNumClutdummy].hingmat!=6 && 
	   nNumClut+1+nNumClutdummy<prot.DOF ;
	 ++nNumClutdummy) {
      //      nNumAtomALL=prot.num_atom-clust[nNumClut].origin_atom_a+1;
      nNumAtomALL+=clust[nNumClut+1+nNumClutdummy].num_atom_clust;
    }
    //    nNumAtomALL+=clust[nNumClut+1+nNumClutdummy].num_atom_clust;
    nNumClutALL=nNumClut/*+1*/+nNumClutdummy;
    if (ABI[nNumClut].hingmat!=6)
      sub_trans_CN_to_A_quaternion(nNumClut,nNumAtomALL,nNumClutALL,/*0*/nNumParent);
    else {
      sub_trans_CN_to_A_quaternion_six(nNumClut,nNumAtomALL,nNumClutALL,/*0*/nNumParent);
    }
    /********************************************************************************/
    /* sub_trans_CN_to_A(nNumClut,prot.num_atom-clust[nNumClut].origin_atom_a+1,0); */
    /********************************************************************************/
  }
	
}


// 局所座標系→実験室系の変換を行う関数
void sub_trans_CN_to_A(int nNumClt, int nNumAtom, int nNumAtomLoca)
{
	int i,j,k,alpha;

	int n=0, angle;

	int nNumClut2;

	int n_delta_dihed;

	int nNumAtomAbsoOrig;

	double sn_delta_dihed;
	double cs_delta_dihed;

	double **Coord/*[MAXA][3]*/;
	double **Coord2/*[MAXA][3]*/;
	double Origin_Coord[3];

	double delta_dihed;

	FILE *outtest;

	Coord=(double **)gcemalloc(sizeof(double *)*prot.num_atom);
	Coord2=(double **)gcemalloc(sizeof(double *)*prot.num_atom);
	for (i=0;i<prot.num_atom;++i) {
	  Coord[i]=(double *)gcemalloc(sizeof(double)*3);
	  Coord2[i]=(double *)gcemalloc(sizeof(double)*3);
	}



	// 現ステップでの二面角の変位
	delta_dihed = clust[nNumClt].now_deltadihedang[0];
//	delta_dihed = -0.1;
	//	delta_dihed = -clust[nNumClt].now_deltadihedang[0];
	if (nNumClt == 1 /*|| nNumClt ==*/ /*3*//*2*//*4*//*5*/)
	{
//		delta_dihed = -clust[nNumClt].now_deltadihedang[0];
	}

///////////////////////////////////////////////////////////
//	if (nNumClt == 3)
//	{
//		delta_dihed += 2.0*PI/10.0;
//	}
//	else
//	{
//		delta_dihed = 0.0;
//	}
///////////////////////////////////////////////////////////
/*	n_delta_dihed = (int)(delta_dihed/(2.0*PI));

	if (delta_dihed > 0.0)
	{
		delta_dihed -= (double)n_delta_dihed*2.0*PI;
	}
	else
	{
		delta_dihed += (double)n_delta_dihed*2.0*PI;
	}

	// 二面角の変位の規格化_1
	if ( 0.0 >= delta_dihed )
	{
		for (;;)
		{
			delta_dihed += 2.0*PI;
			if ( delta_dihed >= 0.0 )
			{
				break;
			}
		}
	}

	// 二面角の変位の規格化_2
	if ( delta_dihed >= 2.0*PI)
	{
		for (;;)
		{
			delta_dihed -= 2.0*PI;
			if ( delta_dihed <= 2.0*PI )
			{
				break;
			}
		}
	}
*/
	sn_delta_dihed = sin(delta_dihed);
	cs_delta_dihed = cos(delta_dihed);

	// 原点と始点の取得
	nNumAtomAbsoOrig = clust[nNumClt].origin_atom_a-1;

	// 原点座標の取得
	for(alpha=0;alpha<3;++alpha)
	{
		Origin_Coord[alpha]=prot.coord[nNumAtomAbsoOrig][alpha];
	}

	// 二面角の変位分だけ局所座標の回転
	for(i = 0; i < nNumAtom ; ++i)
	{
//		Coord[i][0] =   cs_delta_dihed
//		               *clust[nNumClt].xoord_clust/*[0]*/[nNumAtomLoca+i][0]
//		              + sn_delta_dihed
//		               *clust[nNumClt].xoord_clust/*[0]*/[nNumAtomLoca+i][1];
//
//		Coord[i][1] = - sn_delta_dihed
//		               *clust[nNumClt].xoord_clust/*[0]*/[nNumAtomLoca+i][0]
//		              + cs_delta_dihed
//		               *clust[nNumClt].xoord_clust/*[0]*/[nNumAtomLoca+i][1];
//
//		Coord[i][2] = clust[nNumClt].xoord_clust/*[0]*/[nNumAtomLoca+i][2];
		Coord[i][0] =   cs_delta_dihed
		               *clust[nNumClt].xoord_clust/*[0]*/[nNumAtomLoca+i][0]
		              - sn_delta_dihed
		               *clust[nNumClt].xoord_clust/*[0]*/[nNumAtomLoca+i][1];

		Coord[i][1] =   sn_delta_dihed
		               *clust[nNumClt].xoord_clust/*[0]*/[nNumAtomLoca+i][0]
		              + cs_delta_dihed
		               *clust[nNumClt].xoord_clust/*[0]*/[nNumAtomLoca+i][1];

		Coord[i][2] = clust[nNumClt].xoord_clust/*[0]*/[nNumAtomLoca+i][2];
	}

//	// 局所座標系→実験室系の変換_変換行列の乗算
//	for(i = 0; i < nNumAtom ; ++i)
//	{
//		Coord2[i][0] =   clust[nNumClt].trans_A_to_CN[0][0][0]*Coord[i][0]
//		               + clust[nNumClt].trans_A_to_CN[0][1][0]*Coord[i][1]
//		               + clust[nNumClt].trans_A_to_CN[0][2][0]*Coord[i][2];
//
//		Coord2[i][1] =   clust[nNumClt].trans_A_to_CN[0][0][1]*Coord[i][0]
//		               + clust[nNumClt].trans_A_to_CN[0][1][1]*Coord[i][1]
//		               + clust[nNumClt].trans_A_to_CN[0][2][1]*Coord[i][2];
//
//		Coord2[i][2] =   clust[nNumClt].trans_A_to_CN[0][0][2]*Coord[i][0]
//		               + clust[nNumClt].trans_A_to_CN[0][1][2]*Coord[i][1]
//		               + clust[nNumClt].trans_A_to_CN[0][2][2]*Coord[i][2];
//	}

	for(i = 0; i < nNumAtom; ++i)
	{
		for(j = 0; j < 3 ; ++j)
		{
			Coord2[i][j] = 0.0;
		}
	}
	// 局所座標系→実験室系の変換_変換行列の乗算
	for(i = 0; i < nNumAtom; ++i)
	{
		for(j = 0; j < 3 ; ++j)
		{
			for(k = 0; k < 3 ; ++k)
			{
				Coord2[i][j] += clust[nNumClt].trans_A_to_CN[0][k][j]*Coord[i][k];
//				Coord2[i][j] += clust[nNumClt].trans_A_to_CN[0][j][k]*Coord[i][k];
			}
		}
	}

	// 局所座標系→実験室系の変換_原点の移動
	for(i = 0; i < nNumAtom ; ++i)
	{
		for(alpha=0 ;alpha<3; ++alpha)
		{
			prot.coord[nNumAtomAbsoOrig+i][alpha]
			                      = Coord2[i][alpha] + Origin_Coord[alpha];
		}
	}

	// 次のステップの座標での二面角の設定
	clust[nNumClt].dihedang[0] += delta_dihed;

	// 二面角の規格化_1
	if ( 0.0 >= clust[nNumClt].dihedang[0] )
	{
		for (;;)
		{
			clust[nNumClt].dihedang[0] += 2.0*PI;
			if ( clust[nNumClt].dihedang[0] >= 0.0 )
			{
				break;
			}
		}
	}

	// 二面角の規格化_2
	if ( clust[nNumClt].dihedang[0] >= 2.0*PI)
	{
		for (;;)
		{
			clust[nNumClt].dihedang[0] -= 2.0*PI;
			if ( clust[nNumClt].dihedang[0] <= 2.0*PI )
			{
				break;
			}
		}
	}
	
	// 次のステップの座標での局所座標系の計算
	for(nNumClut2=1; nNumClut2<prot.DOF; ++nNumClut2) {
	  trans_A_to_CN(nNumClut2);
	}

//	if ((outtest=fopen("dihedang.out","a")) == NULL)
//	{
//		printf("in\n");
//		exit(1);
//	}
//
//	fprintf(outtest, "%d %e \n",nNumClt, clust[nNumClt].dihedang[0]);
//
//	fclose(outtest);
}

// 局所座標系→実験室系の変換を行う関数
void sub_trans_CN_to_A_quaternion(int nNumClt,int nNumAtomALL, int nNumClutLast, int nNumCltPt) {
  int i,j,k,l;
  int nNumClt2,nNumClut2;
  int nNumAtom;
  int nNumAtomO;
  int nNumAtomP/*,nNumAtomALL*/;
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
  /********************************************************************************************/
  /* nNumAtomALL = clust[nNumClutLast-1].terminal_atom_a[0]-1-clust[nNumClt].origin_atom_a-1; */
  /********************************************************************************************/

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

  /*********************************/
  /* axis_of_rotation_unit[0]=0.0; */
  /* axis_of_rotation_unit[1]=1.0; */
  /* axis_of_rotation_unit[2]=0.0; */
  /* delta_dihed=PI*0.5;	   */
  /*********************************/


  q[0]=cos(delta_dihed*0.5);
  for(i=0;i<3;++i) {
    q[i+1] = axis_of_rotation_unit[i]*sin(delta_dihed*0.5);
  }


  coord_now[0]=0.0;
  for (nNumAtom=0;nNumAtom<nNumAtomALL;++nNumAtom) {
    for (i=0;i<3;++i) {
      coord_now[i+1] = prot.coord[nNumAtom+nNumAtomO][i]-prot.coord[nNumAtomO][i];
    }
    /*********************/
    /* coord_now[0]=0.0; */
    /* coord_now[1]=0.0; */
    /* coord_now[2]=0.0; */
    /* coord_now[3]=1.0; */
    /*********************/
    quaternion_rotation(q,coord_now,coord_rotated);
    //    qua_rot(coord_now,q,coord_rotated);
    for (i=0;i<3;++i) {
      prot.coord[nNumAtom+nNumAtomO][i]= coord_rotated[i+1]+prot.coord[nNumAtomO][i];
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
    

    /**********************************************/
    /* for (i=0;i<3;++i) {			  */
    /*   for (j=0;j<3;++j) {			  */
    /* 	printf("%10.8lf ",dbmat[i][j]);		  */
    /*   }					  */
    /*   printf("\n");				  */
    /* }					  */
    /**********************************************/

  /*******************************/
  /* for (i=0;i<3;++i) {	 */
  /*   for (j=0;j<3;++j) {	 */
  /*     dbmattemp[i][j] = 0.0;	 */
  /*     dbmattemp2[i][j] = 0.0; */
  /*     dbmat2[i][j] = 0.0;	 */
  /*   }			 */
  /* }				 */
  /*******************************/



    /*****************/
    /* printf("\n"); */
    /*****************/


    

  /*******************************************************************************************/
  /* dbmattemp[0][0] = cos(delta_dihed);						     */
  /* dbmattemp[0][1] =-sin(delta_dihed);						     */
  /* dbmattemp[1][0] = sin(delta_dihed);						     */
  /* dbmattemp[1][1] = cos(delta_dihed);						     */
  /* dbmattemp[2][2] = 1.0;								     */
  /*   											     */
  /* for(i = 0; i < 3; ++i) {								     */
  /*   for(j = 0; j < 3 ; ++j) {							     */
  /*     for(k = 0; k < 3 ; ++k) {							     */
  /* 	dbmattemp2[i][j] += clust[nNumClt].trans_A_to_CN[0][k][i]*dbmattemp[k][j];	     */
  /*     }										     */
  /*   }										     */
  /* }											     */
  /*   											     */
  /* for(i = 0; i < 3; ++i) {								     */
  /*   for(j = 0; j < 3 ; ++j) {							     */
  /*     for(k = 0; k < 3 ; ++k) {							     */
  /* 	dbmat2[i][j] += dbmattemp2[i][k]*clust[nNumClt].trans_A_to_CN[0][k][j];		     */
  /*     }										     */
  /*   }										     */
  /* }											     */
  /*******************************************************************************************/
    /* 																	 */
    /* for (nNumAtom=0;nNumAtom<nNumAtomALL;++nNumAtom) {										 */
    /*   																 */
    /*   for (i=0;i<3;++i) {														 */
    /* 	coord_now[i] = prot.coord[nNumAtom+clust[nNumClt].origin_atom_a-1][i]-prot.coord[clust[nNumClt].origin_atom_a-1][i];		 */
    /*   }																 */
    /*   coord_now[3]=0.0;														 */
    /*   for(i = 0; i < 3 ; ++i) 													 */
    /* 	coord_rotated[i] = 0.0;														 */
    /*   																 */
    /*   for(i = 0; i < 3; ++i) {													 */
    /* 	for(j = 0; j < 3 ; ++j) {													 */
    /* 	  coord_rotated[i] += dbmat2[i][j]*coord_now[j];										 */
    /* 	}																 */
    /*   }																 */
    /*   for (i=0;i<3;++i) {														 */
    /* 	prot.coord[nNumAtom+clust[nNumClt].origin_atom_a-1][i]=coord_rotated[i]+prot.coord[clust[nNumClt].origin_atom_a-1][i];		 */
    /*   }																 */
    /* }																 */
    /*************************************************************************************************************************************/
      /****************************************************************************************************************************/

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
    /************************************************************************************/
    /* if (nNumClt2==7) {							        */
    /*   for (i=0;i<3;++i) {							        */
    /* 	for (j=0;j<3;++j) {							        */
    /* 	  printf("%10.8lf ",temp[i][j]);					        */
    /* 	}									        */
    /* 	printf("\n");								        */
    /*   }									        */
    /*   printf("\n");								        */
    /* }									        */
 
 
    for (i=0;i<3;++i) {
      for (j=0;j<3;++j) {
    	clust[nNumClt2].trans_A_to_CN[0][i][j] = temp[i][j];
      }
    }
  }


  /****************************************************/
  /* // 次のステップの座標での局所座標系の計算	      */
  /* 						      */
  /* for(nNumClt2=1; nNumClt2<prot.DOF; ++nNumClt2) { */
  /*   trans_A_to_CN(nNumClt2);			      */
  /* }						      */
  /****************************************************/
  /*   							     */
  /* for (i=0;i<3;++i) {				     */
  /*   for (j=0;j<3;++j) {				     */
  /*     printf("%10.8lf ",clust[7].trans_A_to_CN[0][i][j]); */
  /*   }						     */
  /*   printf("\n");					     */
  /* }							     */
  /* printf("\n");					     */
  /***********************************************************/


    /***********************************************/
    /* for (i=0;i<3;++i) {			   */
    /*   for (j=0;j<3;++j) {			   */
    /* 	printf("%10.8lf ",dbmat2[i][j]);	   */
    /*   }					   */
    /*   printf("\n");				   */
    /* }					   */
    /* printf("\n");				   */
    /***********************************************/

  // 次のステップの座標での二面角の設定
  clust[nNumClt].dihedang[0] += delta_dihed;


    // 二面角の規格化_1
    if ( 0.0 >= clust[nNumClt].dihedang[0] ) {
      for (;;) {
	clust[nNumClt].dihedang[0] += 2.0*PI;
	if ( clust[nNumClt].dihedang[0] >= 0.0) {
	  break;
	}
      }
    }
    
    // 二面角の規格化_2
    if ( clust[nNumClt].dihedang[0] >= 2.0*PI){
      for (;;) {
	clust[nNumClt].dihedang[0] -= 2.0*PI;
	if ( clust[nNumClt].dihedang[0] <= 2.0*PI ) {
	  break;
	}
      }
    }
    
}

void sub_trans_CN_to_A_term(double q[4]) {
  int i,j,k,l,protatomnum,numatomorg;
  double vel[3],acc[3],omg[3],rotmat[3][3];
  double *newcoord,*vec,lenq;
  double dq[4],r[4],roted[4],delta[3];
  double temp[3][3],dbmat[3][3];
  int nNumClut,alpha,alpha2;
  double momentum[7];
  FILE *output,*output2;
  double IV0[6],IV1[6],TIV1[6],tmomentum;


  protatomnum = prot.num_atom;
  if ((newcoord=(double *)calloc(protatomnum*3,sizeof(double))) == NULL) {
    printf("error : alocate \n");
  }
  if ((vec=(double *)calloc(protatomnum*3,sizeof(double))) == NULL) {
    printf("error : alocate \n");
  }

  numatomorg=clust[0].origin_atom_a-1;
  for (i = 0; i < protatomnum; ++i) {
    for (j = 0; j < 3; ++j) {
      vec[i*3+j]=prot.coord[i][j]-prot.coord[numatomorg][j];
    }
  }

  for (i = 0; i < 3; ++i)
    omg[i]=clust[0].sp_velo[i];
  qua_trans_omgtodqua(omg,q,dq);
  for (i=0;i<4;++i)
    q[i]+=deltat*dq[i];

  lenq=0.0;
  for (i=0;i<4;++i)
    lenq+=q[i]*q[i];
  lenq=sqrt(lenq);
  for (i=0;i<4;++i)
    q[i]=q[i]/lenq;

  for (i=0;i<protatomnum;++i) {
    r[0]=0.0;
    for (j=0;j<3;++j)
      r[j+1]=vec[i*3+j];
    quaternion_rotation(q,r,roted);

    for (j=0;j<3;++j)
      newcoord[i*3+j]=prot.coord[numatomorg][j]+roted[j+1];
  }

  for (i = 0; i < 3; ++i) {
    acc[i]=0.0;   	
    vel[i]=0.0;   	
    for (j = 0; j < 3; ++j)
      rotmat[i][j]=0.0;
  }
  qua_set_rotmat(q,rotmat);
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      acc[i]+=rotmat[j][i]*clust[0].sp_acc[j];
      vel[i]+=rotmat[j][i]*clust[0].sp_velo[j];
    }
  }
  
  for(i=0;i<=prot.DOF;++i) {
    for (j=0;j<3;++j) {
      for (k=0;k<3;++k) {
  	temp[j][k] = 0.0;
      }
    }
    for (j=0;j<3;++j) {
      for (k=0;k<3;++k) {
  	for (l=0;l<3;++l) {
  	  temp[j][k] += rotmat[j][k]*clust[i].trans_A_to_CN[0][k][l];
  	}
      }
    }
  
    for (j=0;j<3;++j) {
      for (k=0;k<3;++k) {
	//    	clust[i].trans_A_to_CN[0][j][k] = temp[j][k];
      }
    }
  }

  for (i = 0; i < 3; ++i)
    delta[i]=deltat*vel[i]+0.5*deltat*deltat*acc[i];
  
  for (i=0;i<protatomnum;++i)
    for (j=0;j<3;++j)
      newcoord[i*3+j]+=delta[j];

  /*****************************************/
  /* for (i=0;i<protatomnum;++i)	   */
  /*   for (j=0;j<3;++j)		   */
  /*     prot.coord[i][j]=newcoord[i*3+j]; */
  /*****************************************/
  
  //  output=fopen("test","a");
  //  output2=fopen("momentum","a");
  
  for(i=0; i<prot.num_atom; ++i) {
    for(j=0; j<3; ++j){
      //      fprintf(output,"%12.8lf ",newcoord[i*3+j]);
    }
    //    fprintf(output,"\n ");
  }

  /*****************************************************************/
  /*  for (alpha=0; alpha<6; ++alpha)				   */
  /*   momentum[alpha]=0.0;					   */
  /* for(nNumClut=0; nNumClut<prot.DOF; ++nNumClut){		   */
  /*   for (alpha=0; alpha<6; ++alpha){				   */
  /*     for (alpha2=0; alpha2<6; ++alpha2){			   */
  /* 	momentum[alpha]						   */
  /* 	  +=clust[nNumClut].InertiaMatrix[alpha][alpha2]	   */
  /* 	   *clust[nNumClut].sp_velo[alpha2];			   */
  /*     }							   */
  /*   }							   */
  /* }								   */
  /* 								   */
  /* fprintf(output2,"%d ",nNumStep);				   */
  /* for (alpha=0; alpha<6; ++alpha)				   */
  /*   fprintf(output2,"%12.8lf ",momentum[alpha]);		   */
  /* fprintf(output2,"\n ");					   */
  /*****************************************************************/

  // test for momentum begin
  /***********************************************************/
  /* output2=fopen("momentum","a");			     */
  /* for (i=0; i<7; ++i)				     */
  /*   momentum[i]=0.0;  				     */
  /* mvult(clust[0].InertiaMatrix,clust[0].sp_velo,IV0,6);   */
  /* mvult(clust[1].InertiaMatrix,clust[1].sp_velo,IV1,6);   */
  /* momentum[6]=IV1[2];				     */
  /* 							     */
  /* mvult(clust[1].TransMatrix[0],IV1,TIV1,6);		     */
  /* 							     */
  /* for (i=0; i<6; ++i){				     */
  /*   momentum[i]=IV0[i]+TIV1[i];			     */
  /* }							     */
  /* 							     */
  /* tmomentum=0.0;					     */
  /* for (i = 0; i <7; ++i)				     */
  /*   tmomentum+=momentum[i];				     */
  /* 							     */
  /* 							     */
  /* fprintf(output2,"%d ",nNumStep);			     */
  /* for (i=0; i<7; ++i)				     */
  /*   fprintf(output2,"%12.8lf ",momentum[i]);		     */
  /* fprintf(output2,"sum= %12.8lf ",tmomentum);	     */
  /* fprintf(output2,"\n ");				     */
  /***********************************************************/
  // test for momentum end

  
  //  fprintf(output,"\n\n ");

  free(newcoord);
  free(vec);

  //  fclose(output);
  //  fclose(output2);
}


/***********************************************************************************/
/* void mvult(double m1[6][6], double v[6], double mv[6], int n){		   */
/*   int i,j,k;									   */
/* 										   */
/*   for (i=0;i<n;++i)								   */
/*       mv[i] = 0.0;								   */
/* 										   */
/*   for (i=0;i<n;++i) {							   */
/*     for (j=0;j<n;++j) {							   */
/* 	mv[i] += m1[i][j]*v[j];							   */
/*     }									   */
/*   }										   */
/* }										   */
/* 										   */
/* void mmult(double m1[6][6], double m2[6][6], double m1m2[6][6], int n) {	   */
/*   int i,j,k;									   */
/* 										   */
/*   for (i=0;i<n;++i)								   */
/*     for (j=0;j<n;++j)							   */
/*       m1m2[i][j] = 0.0;							   */
/* 										   */
/*   for (i=0;i<n;++i) {							   */
/*     for (j=0;j<n;++j) {							   */
/*       for (k=0;k<n;++k)	{						   */
/* 	m1m2[i][j] += m1[i][k]*m2[k][j];					   */
/*       }									   */
/*     }									   */
/*   }										   */
/* }										   */
/***********************************************************************************/

// 局所座標系→実験室系の変換を行う関数
void sub_trans_CN_to_A_quaternion_six(int nNumClt,int nNumAtomALL, int nNumClutLast, int nNumCltPt) {
  int i,j,k,l;
  double **coord_now/*[4]*/,**coord_temp2/*[4]*/,**coord_temp/*[4]*/,**coord_rotated/*[3]*/;

  double p[3];
  double ang[3];
  double Rx[3][3],Ry[3][3],Rz[3][3],Rtemp[3][3];
  double R[3][3];
  double T[4][4];

  double origin[3];

  double temp[3][3];

  int nNumAtom;
  int nNumAtomO;

  int nNumClt2;

  nNumAtomO = clust[nNumClt].origin_atom_a-1;
  
  for (i=0;i<3;++i) {
    p[i] = clust[nNumClt].now_deltadihedang_six[i+3];
  }
  for (i=0;i<3;++i) {
    ang[i] = clust[nNumClt].now_deltadihedang_six[i];
  }

  Rx[0][0]=1.0;
  Rx[0][1]=0.0;
  Rx[0][2]=0.0;
  Rx[1][0]=0.0;
  Rx[1][1]=cos(ang[0]);
  Rx[1][2]=-sin(ang[0]);
  Rx[2][0]=0.0;
  Rx[2][1]=sin(ang[0]);
  Rx[2][2]=cos(ang[0]);

  Ry[0][0]=cos(ang[1]);
  Ry[0][1]=0.0;
  Ry[0][2]=-sin(ang[1]);
  Ry[1][0]=0.0;
  Ry[1][1]=1.0;
  Ry[1][2]=0.0;
  Ry[2][0]=sin(ang[1]);
  Ry[2][1]=0.0;
  Ry[2][2]=cos(ang[1]);

  Rz[0][0]=cos(ang[2]);
  Rz[0][1]=-sin(ang[2]);
  Rz[0][2]=0.0;
  Rz[1][0]=sin(ang[2]);
  Rz[1][1]=cos(ang[2]);
  Rz[1][2]=0.0;
  Rz[2][0]=0.0;
  Rz[2][1]=0.0;
  Rz[2][2]=1.0;

  for (i=0;i<3;++i) {
    for (j=0;j<3;++j) {
      Rtemp[i][j]=0.0;
      R[i][j]=0.0;
    }
  }
  
  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	Rtemp[i][j]+=Rx[i][k]*Ry[k][j];

  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	R[i][j]+=Rtemp[i][k]*Rz[k][j];

  for (i=0;i<4;++i)
    for (j=0;j<4;++j)
      T[i][j]=0.0;

  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      T[i][j]=R[i][j];

  for (i=0;i<3;++i)
    T[i][3]=p[i];

  T[3][3]=1.0;

  coord_now=(double **)gcemalloc(sizeof(double*)*nNumAtomALL);
  for (i=0;i<nNumAtomALL;++i)
    coord_now[i]=(double *)gcemalloc(sizeof(double)*4);
  coord_temp2=(double **)gcemalloc(sizeof(double*)*nNumAtomALL);
  for (i=0;i<nNumAtomALL;++i)
    coord_temp2[i]=(double *)gcemalloc(sizeof(double)*4);
  coord_temp=(double **)gcemalloc(sizeof(double*)*nNumAtomALL);
  for (i=0;i<nNumAtomALL;++i)
    coord_temp[i]=(double *)gcemalloc(sizeof(double)*4);
  coord_rotated=(double **)gcemalloc(sizeof(double*)*nNumAtomALL);
  for (i=0;i<nNumAtomALL;++i)
    coord_rotated[i]=(double *)gcemalloc(sizeof(double)*3);
      
  for (nNumAtom=0;nNumAtom<nNumAtomALL;++nNumAtom) {
    coord_now[nNumAtom][3]=1.0;
    for (i=0;i<3;++i) {
      coord_now[nNumAtom][i] = prot.coord[nNumAtom+nNumAtomO][i]-prot.coord[nNumAtomO][i];
    }
    for (i=0;i<3;++i)
      coord_temp[nNumAtom][i]=0.0;
    for (i=0;i<3;++i)
      for (j=0;j<3;++j)
	coord_temp[nNumAtom][i]+=clust[nNumClt].trans_A_to_CN[0][j][i]*coord_now[nNumAtom][j];
    coord_temp[nNumAtom][3]=1.0;

    for (i=0;i<3;++i)
      coord_temp2[nNumAtom][i]=0.0;
    for (i=0;i<4;++i)
      for (j=0;j<4;++j)
	coord_temp2[nNumAtom][i]+=T[i][j]*coord_temp[nNumAtom][j];
    coord_temp2[nNumAtom][3]=1.0;

    for (i=0;i<3;++i)
      coord_rotated[nNumAtom][i]=0.0;
    for (i=0;i<3;++i)
      for (j=0;j<3;++j)
	coord_rotated[nNumAtom][i]+=clust[nNumClt].trans_A_to_CN[0][i][j]*coord_temp2[nNumAtom][j];

  }

  for (i=0;i<3;++i) 
    origin[i]=prot.coord[nNumAtomO][i];
  for (nNumAtom=0;nNumAtom<nNumAtomALL;++nNumAtom) {
    for (i=0;i<3;++i) {
      prot.coord[nNumAtom+nNumAtomO][i]= coord_rotated[nNumAtom][i]+origin[i];
    }
  }

  // 次のステップの座標での局所座標系の計算
  for(nNumClt2=nNumClt;nNumClt2<=nNumClutLast; ++nNumClt2) {
    for (i=0;i<3;++i) {
      for (j=0;j<3;++j) {
	for (k=0;k<3;++k) {
	  temp[i][j] += clust[nNumClt2].trans_A_to_CN[0][i][k]*R[j][k];
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
