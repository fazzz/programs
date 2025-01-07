#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include "clapack.h"

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "force.h"
#include "MD.h"

/**********************************************/
/* void Calc_acc_debug(void);		      */
/* 					      */
/* void mmmult(double *m1, double *m2,int n); */
/* void mvmult(double *m, double *v, int n);  */
/* void mtrans(double *m1, int n);	      */
/* void msetIni(double *m, int n);	      */
/* void msetzero(double *m, int n);	      */
/* int invm(double mat[12][12]);	      */
/**********************************************/

FILE *test;

// 現在の構造でのタンパク質の
// spatial acceleration の計算の補助を行う関数
void sub_calc_sp_acc_cycle(int nNumClut, int nNumCltminusone);

// 現在の構造でのタンパク質の
// spatial acceleration の計算を行う関数
void calc_sp_acc_cycle(int nNumClut){
  int i;
  int nNumClutOfParent;

  nNumClutOfParent = clust[nNumClut].nNumClutOfParent-1;

//	test = fopen("testspac.out", "a");
  // spatial acceleration の初期化を行う
  for(i=0; i<6; ++i){
    clust[nNumClut].sp_acc[i] = 0.0;
  }

//	if (clust[nNumClut].terminal != TERMINAL)
//	{
		// spatial acceleration の計算を行う
  sub_calc_sp_acc_cycle(nNumClut, nNumClutOfParent);
//	}
//	else
//	{
		// spatial acceleration の計算を行う
//		sub_sub_calc_sp_acc_cycle(nNumClut);
//	}
//	fclose(test);
		
  if (TermMoveMode==ON){
    if (nNumClut==0){
      for (i=0;i<6;++i) {
	clust[0].sp_acc[i]=acc_Term[i];
      }
    }
    if (nNumClut==0 && TermMoveMode2 == 12){
      for (i=0;i<6;++i) {
    	clust[0].sp_acc[i]=acc_Term2[i];
      }
    }
  }

  /******************************/
  /* if (nNumClut==1)	      */
  /*   Calc_acc_debug();	      */
  /******************************/

}

// 現在の構造でのタンパク質の
// spatial acceleration の計算の補助を行う関数
void sub_calc_sp_acc_cycle(int nNumClut, int nNumCltminusone) {
  int i,j;

  int alpha,alpha2;

  int nNumFreedom;

  double DotOmegaOfThisCoord[3];
  double DotOmegaOfAbsoCoord[3];
  double OmegaOmegaL[3];

  double q[3];

  FILE *outtest;
  FILE *debug;

  // 剛体のspatial accerlatontの初期化を行う
  for(i=0; i<6; ++i)  {
    clust[nNumClut].sp_acc[i] = 0.0;
  }

  // 剛体の予測加速度の初期化を行う
  for(i=0; i<6; ++i) {
    clust[nNumClut].predict_alpha[i] = 0.0;
  }

  // 変換行列の転置行列を乗じる
  for (i=0;i<6;++i)  {
    for (j=0;j<6;++j)  {
      clust[nNumClut].predict_alpha[i]
	+=   clust[nNumClut].TransMatrix[0][j][i]
	* clust[nNumCltminusone].sp_acc[j];
    }
  }

	// Coriolis 加速度部分の初期化
//	for (alpha=0;alpha<3;++alpha)
//	{
//		OmegaOmegaL[alpha] = 0.0;
//	}

  // 剛体の2面角加速度の初期化を行う
  clust[nNumClut].dddihedang[0] = 0.0;


////////////////////////////////////////////////////////////////////////
	// 剛体の2面角加速度の計算を行う_1
  for(i=0; i<6; ++i) {
    clust[nNumClut].dddihedang[0]/*rad/s^2*/
      -=   ABI[nNumClut].Kalman_Gain[i]
      *clust[nNumClut].predict_alpha[i]/*rad/s^2*/;
  }

  // 剛体の2面角加速度の計算を行う_2
  clust[nNumClut].dddihedang[0]/*rad/s^2*/+= zzz[nNumClut].nyu/*rad/s^2*/;

//	if (MODE == NVT)
//	{
//		s_NVT = s_NVT_predict[0];
//		dot_s_NVT = s_NVT_predict[1]/deltat;
//		clust[nNumClut].dddihedang[0] -= dot_s_NVT/s_NVT/*xi_dummy**//*xi_dummy**/*clust[nNumClut].ddihedang[0];
//	}

//	pick_dihed2(nNumClut);

//		if ((outtest=fopen("dddihedang.out","a")) == NULL)
//		{
//			printf("in\n");
//			exit(1);
//		}
//		for(nNumClut=1; nNumClut<prot.DOF; ++nNumClut)
//		{
//			fprintf(outtest, "%d %e %e %e \n",nNumClut,clust[nNumClut].dihedang[0],clust[nNumClut].ddihedang[0],clust[nNumClut].dddihedang[0]);
//		}
//		fclose(outtest);

//	fprintf(test,"dddihedang[%d] = %e \n", nNumClut, clust[nNumClut].dddihedang[0]/**1.0e20*/);
////////////////////////////////////////////////////////////////////////

	// spatial acceleration の計算を行う_1
  for(i=0; i<6; ++i)	{
    clust[nNumClut].sp_acc[i]/*rad/s^2*/
      +=   clust[nNumClut].predict_alpha[i]/*rad/s^2*/
      + clust[nNumClut].Coriolis_acc[i];/*rad/s^2*/
  }

	// クラスタ座標系での角加速度と
	// 実験室座標系での角加速度の初期化
//	for (alpha=0;alpha<3;++alpha)
//	{
//		DotOmegaOfThisCoord[alpha] = 0.0;
//		DotOmegaOfAbsoCoord[alpha] = 0.0;
//	}

	// このクラスタの自由度の取得
  nNumFreedom = ABI[nNumClut].hingmat-1;

	// spatial acceleration の計算を行う_2
  clust[nNumClut].sp_acc[nNumFreedom] += clust[nNumClut].dddihedang[0];


//	fprintf(test,"sp_acc[%d] = ", nNumClut);
//	for (i=0;i<6;++i)
//	{
//		fprintf(test,"%e \n", clust[nNumClut].sp_acc[i]/**1.0e20*/);
//	}
//	fprintf(test, "\n");
	// クラスタ座標系での角加速度の取得
//	DotOmegaOfThisCoord[nNumFreedom] = clust[nNumClut].dddihedang[0];

	// 実験室座標系での角加速度の計算
//	for (alpha=0;alpha<3;++alpha)
//	{
//		for (alpha2=0;alpha2<3;++alpha2)
//		{
//			DotOmegaOfAbsoCoord[alpha]
//				+= clust[nNumClut].trans_A_to_CN[0][alpha2][alpha]
//				  *DotOmegaOfThisCoord[alpha2];
//		}
//	}

	 // 実験室座標系での角速度
//	DotOmegaOfAbsoCoord[0]
//		=   clust[nNumClut].trans_A_to_CN[0][1][1]
//		   *clust[nNumClut].trans_A_to_CN[0][2][0]
//		   *clust[nNumClut].ddihedang[0]
//		  - clust[nNumClut].trans_A_to_CN[0][1][0]
//		   *clust[nNumClut].trans_A_to_CN[0][2][1]
//		   *clust[nNumClut].ddihedang[0];

//	DotOmegaOfAbsoCoord[1]
//		=   clust[nNumClut].trans_A_to_CN[0][0][1]
//		   *clust[nNumClut].trans_A_to_CN[0][2][0]
//		   *clust[nNumClut].ddihedang[0]
//		  - clust[nNumClut].trans_A_to_CN[0][0][0]
//		   *clust[nNumClut].trans_A_to_CN[0][2][1]
//		   *clust[nNumClut].ddihedang[0];

//	DotOmegaOfAbsoCoord[2]
//		=   clust[nNumClut].trans_A_to_CN[0][0][0]
//		   *clust[nNumClut].trans_A_to_CN[0][1][1]
//		   *clust[nNumClut].ddihedang[0]
//		  - clust[nNumClut].trans_A_to_CN[0][0][1]
//		   *clust[nNumClut].trans_A_to_CN[0][1][0]
//		   *clust[nNumClut].ddihedang[0];

	// クラスタ間のベクトルの取得
//	for (alpha=0;alpha<3;++alpha)
//	{
//		q[alpha] = clust[nNumClut].TransMatrix[0][alpha+3][alpha2];
//	}

	// クラスタ間のベクトルの取得
//	q[0] = clust[nNumClut].TransMatrix[0][1][3];
//	q[1] = clust[nNumClut].TransMatrix[0][0][5];
//	q[2] = clust[nNumClut].TransMatrix[0][2][4];

	// Coriolis 加速度部分の計算_1
//	for (alpha=0;alpha<3;++alpha)
//	{
//		for (alpha2=0;alpha2<3;++alpha2)
//		{
//			OmegaOmegaL[alpha]
//			   += clust[nNumClut].Coriolis_acc_Mat[alpha][alpha2]*q[alpha2];
//		}
//	}

	// Coriolis 加速度部分の計算_2
//	for (alpha=0;alpha<3;++alpha)
//	{
//		clust[nNumClut].sp_acc[alpha+3]/*rad/s^2*/
//		                 += OmegaOmegaL[alpha];/*rad/s^2*/
//	}

	// spatial acceleration の計算を行う_2
//	for (alpha=0;alpha<3;++alpha)
//	{
//		clust[nNumClut].sp_acc[alpha]/*rad/s^2*/ 
//		                 += DotOmegaOfAbsoCoord[alpha];/*rad/s^2*/
//	}
}

// 現在の構造でのタンパク質の
// spatial acceleration の計算の補助を行う関数
//void sub_sub_calc_sp_acc_cycle(int nNumClut)
//{
//	int i,j;
//
//	int alpha,alpha2;
//
//	int nNumFreedom;
//
//	double DotOmegaOfThisCoord[3];
//	double DotOmegaOfAbsoCoord[3];
//	double OmegaOmegaL[3];
//
//	double TransMatrix[6][6];
//
//	double q[3];
//
//	// 剛体のspatial accerlatontの初期化を行う
//	for(i=0; i<6; ++i)
//	{
//		clust[nNumClut].sp_acc[i] = 0.0;
//	}
//
//	// Coriolis 加速度部分の初期化
//	for (alpha=0;alpha<3;++alpha)
//	{
//		OmegaOmegaL[alpha] = 0.0;
//	}
//
//	// 剛体の2面角加速度の初期化を行う
//	clust[nNumClut].dddihedang[0] = 0.0;
//
//////////////////////////////////////////////////////////////////////////
//
//	// 剛体の2面角加速度の計算を行う_1
//	for(i=0; i<6; ++i)
//	{
//		clust[nNumClut].dddihedang[0]/*rad/s^2*/
//						-=   ABI[nNumClut].Kalman_Gain_Transpose[i]
//					    	*clust[nNumClut].predict_alpha[i]/*rad/s^2*/;
//	}
//
//	// 剛体の2面角加速度の計算を行う_2
//	clust[nNumClut].dddihedang[0]/*rad/s^2*/ 
//						+= zzz[nNumClut].nyu/*rad/s^2*/;
//////////////////////////////////////////////////////////////////////////
//
//	// spatial acceleration の計算を行う_1
//	for(i=0; i<6; ++i)
//	{
//		clust[nNumClut].sp_acc[i]/*rad/s^2*/
//		                 +=   clust[nNumClut].predict_alpha[i]/*rad/s^2*/;
//// 			 		        + clust[nNumClut].Coriolis_acc[i];/*rad/s^2*/
//	}
//
//	// クラスタ座標系での角加速度と
//	// 実験室座標系での角加速度の初期化
//	for (alpha=0;alpha<3;++alpha)
//	{
//		DotOmegaOfThisCoord[alpha] = 0.0;
//		DotOmegaOfAbsoCoord[alpha] = 0.0;
//	}
//
//	// このクラスタの自由度の取得
//	nNumFreedom = ABI[nNumClut].hingmat-1;
//
//	// クラスタ座標系での角加速度の取得
//	DotOmegaOfThisCoord[nNumFreedom] = clust[nNumClut].dddihedang[0];
//
//	// 実験室座標系での角加速度の計算
//	for (alpha=0;alpha<3;++alpha)
//	{
//		for (alpha2=0;alpha2<3;++alpha2)
//		{
//			DotOmegaOfAbsoCoord[alpha]
//				+= clust[nNumClut].trans_A_to_CN[0][alpha2][alpha]
//				  *DotOmegaOfThisCoord[alpha2];
//		}
//	}
//
//	for (i=0;i<6;++i)
//	{
//		for (j=0;j<6;++j)
//		{
//			TransMatrix[i][j] = 0.0;
//		}
//	}
//
//	for (i=0;i<6;++i)
//	{
//		for (j=0;j<6;++j)
//		{
//			TransMatrix[i][j] = -1*clust[nNumClut].TransMatrix[0][i][j];
//		}
//	}
//
//	for (i=0;i<6;++i)
//	{
//		TransMatrix[i][i] += 2;
//	}
//
//	// クラスタ間のベクトルの取得
//	for (alpha=0;alpha<3;++alpha)
//	{
//		q[alpha] = TransMatrix[alpha+3][alpha2];
//	}
//
//	// クラスタ間のベクトルの取得
//	q[0] = TransMatrix[1][3];
//	q[1] = TransMatrix[0][5];
//	q[2] = TransMatrix[2][4];
//
//	// Coriolis 加速度部分の計算_1
//	for (alpha=0;alpha<3;++alpha)
//	{
//		for (alpha2=0;alpha2<3;++alpha2)
//		{
//			OmegaOmegaL[alpha]
//			   += clust[nNumClut].Coriolis_acc_Mat[alpha][alpha2]*q[alpha2];
//		}
//	}
//
//	// Coriolis 加速度部分の計算_2
//	for (alpha=0;alpha<3;++alpha)
//	{
////		clust[nNumClut].sp_acc[alpha+3]/*rad/s^2*/
////		                 += OmegaOmegaL[alpha];/*rad/s^2*/
//	}
//
//	// spatial acceleration の計算を行う_2
//	for (alpha=0;alpha<3;++alpha)
//	{
//		clust[nNumClut].sp_acc[alpha]/*rad/s^2*/ 
//		                 += DotOmegaOfAbsoCoord[alpha];/*rad/s^2*/
//	}
//}

/*********************************************************************/
/* void Calc_acc_debug(void){					     */
/*   int i,j;							     */
/*   double *Mat;						     */
/*   double *MatPsi,*MatInt,*MatHing;				     */
/*   double *force,*biasf,*T,*theta;				     */
/*   double MatHingPsiIintPsiTHingT2[12][12];			     */
/* 								     */
/*   Mat=(double *)malloc(sizeof(double)*12*12);		     */
/*   MatPsi=(double *)malloc(sizeof(double)*12*12);		     */
/*   MatInt=(double *)malloc(sizeof(double)*12*12);		     */
/*   MatHing=(double *)malloc(sizeof(double)*12*12);		     */
/*   force=(double *)malloc(sizeof(double)*12);			     */
/*   biasf=(double *)malloc(sizeof(double)*12);			     */
/*   T=(double *)malloc(sizeof(double)*12);			     */
/*   theta=(double *)malloc(sizeof(double)*12);			     */
/* 								     */
/*   msetIni(MatPsi,12);					     */
/*   msetzero(MatInt,12); 					     */
/*   msetIni(MatHing,12);					     */
/*   								     */
/*   for (i=0;i<6;++i) {					     */
/*     for (j=0;j<6;++j) {					     */
/*       MatPsi[(i+6)*12+j]=clust[1].TransMatrix[0][i][j];	     */
/*       MatInt[i*12+j]=clust[0].InertiaMatrix[i][j];		     */
/*       MatInt[(i+6)*12+(j+6)]=clust[1].InertiaMatrix[i][j];	     */
/*     }							     */
/*   }								     */
/*   MatHing[6*12+6]=0.0;					     */
/*   MatHing[7*12+7]=0.0;					     */
/*   MatHing[8*12+8]=1.0;					     */
/*   MatHing[9*12+9]=0.0;					     */
/*   MatHing[10*12+10]=0.0;					     */
/*   MatHing[11*12+11]=0.0;					     */
/* 								     */
/*   force[8]=clust[1].f_c.f_dihed;				     */
/*   for (i=0;i<6;++i) {					     */
/*     biasf[i]=clust[0].Coriolis_b[i];				     */
/*     biasf[6+i]=clust[0].Coriolis_b[i];			     */
/*   }								     */
/*   mvmult(MatHing,force,12);					     */
/*   mvmult(MatPsi,biasf,12);					     */
/*   mvmult(MatHing,biasf,12);					     */
/*   for (i=0;i<12;++i) {					     */
/*     T[i]=force[i]-biasf[i];					     */
/*   }								     */
/*  								     */
/* 								     */
/*   for (i=0;i<12;++i)						     */
/*     for (j=0;j<12;++j)					     */
/*       Mat[i*12+j]=MatInt[i*12+j];				     */
/*   mmmult(MatPsi,Mat,12);					     */
/*   mtrans(MatPsi,12);						     */
/*   mmmult(Mat,MatPsi,12);					     */
/*   for (i=0;i<12;++i)						     */
/*     for (j=0;j<12;++j)					     */
/*       Mat[i*12+j]=MatPsi[i*12+j];				     */
/*   mmmult(MatHing,Mat,12);					     */
/*   mmmult(Mat,MatHing,12);					     */
/*   for (i=0;i<12;++i)						     */
/*     for (j=0;j<12;++j)					     */
/*       Mat[i*12+j]=MatHing[i*12+j];				     */
/*   for (i=0;i<12;++i)						     */
/*     for (j=0;j<12;++j)					     */
/*       MatHingPsiIintPsiTHingT2[i][j]=Mat[i*12+j];		     */
/*   invm(MatHingPsiIintPsiTHingT2);				     */
/*   for (i=0;i<12;++i)						     */
/*     for (j=0;j<12;++j)					     */
/*       Mat[i*12+j]=MatHingPsiIintPsiTHingT2[i][j];		     */
/* 								     */
/* 								     */
/*  mvmult(Mat,T,12);						     */
/*   for (i=0;i<12;++i) {					     */
/*     theta[i]=T[i];						     */
/*   }								     */
/* 								     */
/*   free(MatPsi);						     */
/*   free(MatInt);						     */
/*   free(MatHing);						     */
/*   free(Mat);							     */
/*   free(force);						     */
/*   free(biasf);						     */
/*   free(T);							     */
/*   free(theta);						     */
/* 								     */
/* 								     */
/* }								     */
/* 								     */
/* void mmmult(double *m1, double *m2, int n){			     */
/*   int i,j,k;							     */
/*   double *m1m2;						     */
/* 								     */
/*   m1m2=(double *)malloc(sizeof(double)*n*n);			     */
/*   								     */
/*   for (i=0;i<n;++i)						     */
/*     for (j=0;j<n;++j)					     */
/*       m1m2[i*n+j] = 0.0;					     */
/* 								     */
/*   for (i=0;i<n;++i)						     */
/*     for (j=0;j<n;++j)					     */
/*       for (k=0;k<n;++k)					     */
/* 	m1m2[i*n+j] += m1[i*n+k]*m2[k*n+j];			     */
/* 								     */
/*   for (i=0;i<n;++i)						     */
/*     for (j=0;j<n;++j)					     */
/*       m2[i*n+j] = m1m2[i*n+j];				     */
/* 								     */
/* 								     */
/*   free(m1m2);						     */
/* }								     */
/* 								     */
/* void mvmult(double *m, double *v, int n){			     */
/*   int i,j,k;							     */
/*   double *mv;						     */
/* 								     */
/*   mv=(double *)malloc(sizeof(double)*n);			     */
/*   								     */
/*   for (i=0;i<n;++i)						     */
/*       mv[i] = 0.0;						     */
/* 								     */
/*   for (i=0;i<n;++i)						     */
/*     for (j=0;j<n;++j)					     */
/* 	mv[i] += m[i*n+j]*v[j];					     */
/* 								     */
/*   for (i=0;i<n;++i)						     */
/*     v[i] = mv[i];						     */
/* 								     */
/* 								     */
/*   free(mv);							     */
/* 								     */
/* }								     */
/* 								     */
/* void mtrans(double *m1,int n){				     */
/*   int i,j;							     */
/*   double *m2;						     */
/* 								     */
/*   m2=(double *)malloc(sizeof(double)*n*n);			     */
/* 								     */
/*   for (i=0;i<n;++i)						     */
/*     for (j=0;j<n;++j)					     */
/*       m2[i*n+j] = m1[j*n+i];					     */
/* 								     */
/*   for (i=0;i<n;++i)						     */
/*     for (j=0;j<n;++j)					     */
/*       m1[i*n+j] = m2[i*n+j];					     */
/* 								     */
/*   free(m2);							     */
/* }								     */
/* 								     */
/* void msetIni(double *m, int n) {				     */
/*   int i,j;							     */
/* 								     */
/* 								     */
/*   for (i=0;i<n;++i)						     */
/*     for (j=0;j<n;++j)					     */
/*       m[i*n+j] = 0.0;					     */
/* 								     */
/*   for (i=0;i<n;++i)						     */
/*       m[i*n+i] = 1.0;					     */
/* 								     */
/* }								     */
/* 								     */
/* void msetzero(double *m, int n) {				     */
/*   int i,j;							     */
/* 								     */
/* 								     */
/*   for (i=0;i<n;++i)						     */
/*     for (j=0;j<n;++j)					     */
/*       m[i*n+j] = 0.0;					     */
/* 								     */
/* }								     */
/* 								     */
/* int invm(double mat[12][12]) {				     */
/* 								     */
/*   int i,j,k;							     */
/*   double A[144];						     */
/*   double inv[12][12];					     */
/*   double test[12][12];					     */
/*   static long int m=12,n=12,lda=12,info,piv[12],lwork=12;	     */
/*   static double work[12];					     */
/* 								     */
/*   k=0;							     */
/*   for(i=0;i<12;++i) {					     */
/*     for(j=0;j<12;++j) {					     */
/*       A[k]=mat[j][i];					     */
/*       ++k;							     */
/*     }							     */
/*   }								     */
/* 								     */
/*   dgetrf_(&m,&n,A,&lda,piv,&info);				     */
/* 								     */
/*   dgetri_(&n,A,&lda,piv,work,&lwork,&info);			     */
/* 								     */
/*   k=0;							     */
/*   for(i=0;i<12;++i) {					     */
/*     for(j=0;j<12;++j) {					     */
/*       inv[j][i]=A[k];					     */
/*       ++k;							     */
/*     }							     */
/*   }								     */
/* 								     */
/*   for(i=0;i<12;++i)						     */
/*     for(j=0;j<12;++j)					     */
/*       test[i][j] = 0.0;					     */
/* 								     */
/*   for(i=0;i<12;++i)						     */
/*     for(j=0;j<12;++j)					     */
/*       for(k=0;k<12;++k)					     */
/* 	test[i][j]+=mat[i][k]*inv[k][j];			     */
/* 								     */
/*   for(i=0;i<12;++i)						     */
/*     for(j=0;j<12;++j)					     */
/*       mat[i][j] = inv[i][j];					     */
/* 								     */
/*   return 0;							     */
/* 								     */
/* }								     */
/*********************************************************************/
