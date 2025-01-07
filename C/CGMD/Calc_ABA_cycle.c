#include <stdio.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "MD.h"
#include "BD.h"

double MatHing[6];
double MatHingTranspose[6];

//FILE *outtest;

// Articulated-Body Inertia 、
// Kalman Gain の計算を行う関数
void calc_ABA_cycle(int nNumClut)
{
	int nNumTerminal;
	int nNumClutTarget;
	int nNumClutOfChild[10];
	int child;
	nNumClutTarget = IndexOfABICycle[nNumClut]- 1;

//	outtest = fopen("testACI.out", "a");
	// 末端の ABI 、KG の計算
//	if (nNumClutTarget==prot.DOF-1)
//	if (clust[nNumClutTarget].terminal == 0)
//	{
//		sub_calc_initial_ABI(nNumClutTarget);
//	}
//	else
//	{
	if (nNumClutTarget!=-1) { // 0811
		for (child=0;child<clust[nNumClutTarget].num_branch;++child)
		{
			nNumClutOfChild[child]
			      = clust[nNumClutTarget].nNumClutOfChild[child]-1;
		}
	}                        // 0811
		sub_calc_ABI_cycle(nNumClutTarget ,
			               nNumClutOfChild,
						   clust[nNumClutTarget].num_branch);
//	}
//	fclose(outtest);
}

// 末端の ABI 、KG の計算を行う関数
void sub_calc_initial_ABI(int nNumClut)
{
	int i,j;

	if (nNumClut!=-1) { // 0811
//	fprintf(outtest,"CorrectABI[prot.DOF-1]\n");
	for(i=0;i<6;++i)
	{
		for(j=0;j<6;++j)
		{
			ABI[nNumClut].CorrectABI[i][j]
				= clust[nNumClut].InertiaMatrix[i][j]/*kg*m^2/kg/kgm*/;
//			fprintf(outtest,"%e ",ABI[prot.DOF-1].CorrectABI[i][j]/**1.0e40*/);
//			ABI[prot.DOF-1].PredictABI[i][j]
//				= clust[prot.DOF-1].InertiaMatrix[i][j]/*kg*m^2/kg/kgm*/;
		}
//			fprintf(outtest,"\n");
	}
	} // 0811
//	fprintf(outtest,"\n\n");
	CalcD(nNumClut);
	calcKalmanGain(nNumClut);
	calcPredictABI(nNumClut);
}

// ABI の計算を行う関数
void sub_calc_ABI_cycle(int nNumCluto, // 現在の剛体のインデックス
	                    int nNumClutplusone[10], // 現在の剛体の親剛体のインデックス
                        int num_branch)
{
	// 修正子 ABI の計算を行う
	CalcCorrectABI(nNumCluto, nNumClutplusone, num_branch);
	// Joint Inertia の計算を行う
	CalcD(nNumCluto);
	// KG の計算を行う
	calcKalmanGain(nNumCluto);
	// 予測子 ABI の計算を行う
	calcPredictABI(nNumCluto);
}

// ABI の計算を行う関数
void sub_sub_calc_ABI_cycle(int nNumCluto, int nNumClutplusone)
{
	// 修正子 ABI の計算を行う
	sub_CalcCorrectABI(nNumCluto, nNumClutplusone);
	// Joint Inertia の計算を行う
	CalcD(nNumCluto);
	// KG の計算を行う
	calcKalmanGain(nNumCluto);
	// 予測子 ABI の計算を行う
	calcPredictABI(nNumCluto);
}

// 予測子 ABI の計算を行う関数///////////////////////////////////////
int calcPredictABI(int nNumCluto)
{
	int i,j,k;

	int num_dof;

	double tau[6][6];
	double hingmat[6];

	for (i=0;i<6;++i)
	{
		hingmat[i] = 0.0;
	}

	if (nNumCluto!=-1) { // 0811
	// 予測子 ABI の初期化
	for (i=0;i<6;++i)
	{
		for (j=0;j<6;++j)
		{
			ABI[nNumCluto].PredictABI[i][j] = 0.0;
		}
	}
	}
	num_dof=/*ABI[nNumCluto].hingmat-1*/2;

//	hingmat[num_dof] = 1.0;
//
//	for (i=0;i<6;++i)
//	{
//		for (j=0;j<6;++j)
//		{
//			tau[i][j]
//			   =  ident[i][j]
//			       - ( ABI[nNumCluto].Kalman_Gain[i]*hingmat[j]);
//		}
//	}
//
////	fprintf(outtest,"PredictABI %d = \n", nNumCluto);
//	for(i=0;i<6;++i)
//	{
//		for(j=0;j<6;++j)
//		{
//			for(k=0;k<6;++k)
//			{
//				ABI[nNumCluto].PredictABI[i][j]
//				  += tau[i][k]*ABI[nNumCluto].CorrectABI[k][j];
//			}
////			fprintf(outtest,"%e ",ABI[nNumCluto].PredictABI[i][j]/**1.0e40*/);
//		}
////		fprintf(outtest,"\n");
//	}
//	fprintf(outtest,"\n\n");

	if (nNumCluto!=-1) { // 0811
	for (i=0;i<6;++i)
	{
		for (j=0;j<6;++j)
		{
			ABI[nNumCluto].PredictABI[i][j]
			=ABI[nNumCluto].CorrectABI[i][j] - ABI[nNumCluto].CorrectABI[i][2]
											  *ABI[nNumCluto].CorrectABI[2][j]
											  /ABI[nNumCluto].CorrectABI[2][2];
		}
	}
	}
	// BD の場合
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/*	if (DYNMODE == BD)
	{
		for(i=0;i<6;++i)
		{
			for(j=0;j<6;++j)
			{
				ABI[nNumCluto].PredictABI[i][j]
				  = clust[nNumCluto].friction_tensor_tra
				   *ABI[nNumCluto].PredictABI[i][j];
			}
		}
	}*/
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//	// tau の計算
//	for (i=0;i<6;++i)
//	{
//		for (j=0;j<6;++j)
//		{
//			tau[i][j] =   ident[i][j]
//                        - ABI[nNumCluto].Kalman_Gain[i]
//                         *MatHing[j];
//		}
//	}

	// 予測子 ABI の計算
//	for(i=0;i<6;++i)
//	{
//		for(j=0;j<6;++j)
//		{
//			for(k=0;k<6;++k)
//			{
//				ABI[nNumCluto].PredictABI[i][j]
//						 += tau[i][k]*ABI[nNumCluto].CorrectABI[k][j];
//			}
//		}
//	}
//	if (nNumCluto == prot.DOF-1)
//	{
//		for (i=0;i<6;++i)
//		{
//			for (j=0;j<6;++j)
//			{
//				ABI[nNumCluto].PredictABI[i][j] = 0.0;
//			}
//		}
//	}

}
/////////////////////////////////////////////////////////////////////

// KG の計算を行う関数
void calcKalmanGain(int nNumCluto){
  int i,j;
  int num_dof;
  // KG の初期化
  if (nNumCluto!=-1) { // 0811
  for(i=0;i<6;++i){
    ABI[nNumCluto].Kalman_Gain[i] = 0.0;
  }
  num_dof=/*ABI[nNumCluto].hingmat-1*/2;
  // KG の計算を行う
  for(i=0; i<6; ++i) {
    ABI[nNumCluto].Kalman_Gain[i]= ABI[nNumCluto].CorrectABI[i][num_dof]/ABI[nNumCluto].ABJointI;
  }
  }
}

// Joint Inertia の計算を行う関数
void CalcD(int nNumCluto){
  int num_dof;
  if (nNumCluto!=-1) { // 0811
  num_dof=/*ABI[nNumCluto].hingmat-1*/2;
  ABI[nNumCluto].ABJointI = ABI[nNumCluto].CorrectABI[num_dof][num_dof];
  }
}

// 修正子 ABI の計算を行う関数_1
void CalcCorrectABI(int nNumCluto, 
		    int nNumClutplusone[10], 
		    int num_branch){
  int i,j,k,l;
  int child,nNumNow;
  double mat[6][6];

  // 修正子 ABI の初期化
  if (nNumCluto!=-1) { // 0811
  for(i=0;i<6;++i){
    for(j=0;j<6;++j){
      mat[i][j] = 0.0;
      ABI[nNumCluto].CorrectABI[i][j] = 0.0;
    }
  }
  }

  if (nNumCluto==-1) num_branch=0; // 08_11

  for (child=0;child<num_branch;++child){
    nNumNow = nNumClutplusone[child];
    if (nNumCluto!=-1 && nNumNow != -1) { // 0811
    if (nNumNow != -2){
///////////////////////////////////////////////////////////////////////////////
      for(i=0;i<6;++i){
	for(j=0;j<6;++j){
	  for(k=0;k<6;++k){
	    for(l=0;l<6;++l){
	      //						nNumNow = nNumClutplusone[child];
	      // 変換行列、予測子 ABI 、変換行列の転置行列
	      ABI[nNumCluto].CorrectABI[i][j]
		+= clust[nNumNow/*nNumCluto*/].TransMatrix[0][/*k*/i][/*i*/k]
		*ABI[nNumNow].PredictABI[k][l]
		*clust[nNumNow/*nNumCluto*/].TransMatrix[0][/*l*/j][/*j*/l];
	    }
	  }
	}
      }
///////////////////////////////////////////////////////////////////////////////
    }
    }
//		nNumNow = nNumClutplusone[child];
//		for(i=0;i<6;++i)
//		{
//			for(j=0;j<6;++j)
//			{
//				for(k=0;k<6;++k)
//				{
//					// 変換行列、予測子 ABI 、変換行列の転置行列
//					mat[i][j]
//						+= clust[nNumNow].TransMatrix[0][i][k]
//						  *ABI[nNumNow].PredictABI[k][j];
//				}
//			}
//		}
//
//		for(i=0;i<6;++i)
//		{
//			for(j=0;j<6;++j)
//			{
//				for(k=0;k<6;++k)
//				{
//					// 変換行列、予測子 ABI 、変換行列の転置行列
//					ABI[nNumCluto].CorrectABI[i][j]
//						+= mat[i][k]
//						  *clust[nNumNow].TransMatrix[0][j][k];
//				}
//			}
//		}
///////////////////////////////////////////////////////////////////////////////
	}

//	for(i=0;i<6;++i)
//	{
//		for(j=0;j<6;++j)
//		{
//			for(k=0;k<6;++k)
//			{
//				for(l=0;l<6;++l)
//				{
//					// 変換行列、予測子 ABI 、変換行列の転置行列
//					ABI[nNumCluto].CorrectABI[i][j]
//						+= clust[nNumCluto].TransMatrix[0][i][k]
//						  *ABI[nNumClutplusone].PredictABI[k][l]
//					      *clust[nNumCluto].TransMatrix[0][j][l];
//				}
//			}
//		}
//	}

//	fprintf(outtest,"CorrectABI %d = ", nNumCluto);
  if (nNumNow != -1) { // 0811
  for(i=0;i<6;++i) {
    for(j=0;j<6;++j){
      // 慣性行列を足す
      ABI[nNumCluto].CorrectABI[i][j] 
	+= clust[nNumCluto].InertiaMatrix[i][j]
	/*kg*m^2/kg/kgm*/;
      //	fprintf(outtest,"%e ",ABI[nNumCluto].CorrectABI[i][j]/**1.0e40*/);
    }
    //	fprintf(outtest,"\n");
  }
  }
//	fprintf(outtest,"\n");
	// BD の場合
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/*	else if(DYNMODE == BD)
	{
		for(i=0;i<6;++i)
		{
			for(j=0;j<6;++j)
			{
				// 慣性行列を足す
				ABI[nNumCluto].CorrectABI[i][j] 
				      += clust[nNumCluto].friction_tensor_tra
						*clust[nNumCluto].InertiaMatrix[i][j]
				      /*kg*m^2/kg/kgm*//*;
			}
		}
	}*/
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
}

// 修正子 ABI の計算を行う関数_2
void sub_CalcCorrectABI(int nNumCluto, int nNumClutplusone)
{
	int i,j,k,l;

	double TransMatrix[6][6];

	if (nNumCluto != -1) { // 0811 

	// 修正子 Bias force の初期化
	for(i=0;i<6;++i)
	{
		zzz[nNumCluto].Correctzzz[i] = 0.0;
	}

	for (i=0;i<6;++i)
	{
		for (j=0;j<6;++j)
		{
			TransMatrix[i][j] = 0.0;
		}
	}

	for (i=0;i<6;++i)
	{
		for (j=0;j<6;++j)
		{
			TransMatrix[i][j] = -1*clust[nNumCluto].TransMatrix[0][i][j];
		}
	}

	for (i=0;i<6;++i)
	{
		TransMatrix[i][i] += 2;
	}

	// 修正子 ABI の初期化
	for(i=0;i<6;++i)
	{
		for(j=0;j<6;++j)
		{
			ABI[nNumCluto].CorrectABI[i][j] = 0.0;
		}
	}

	for(i=0;i<6;++i)
	{
		for(j=0;j<6;++j)
		{
			for(k=0;k<6;++k)
			{
				for(l=0;l<6;++l)
				{
//					// 変換行列、予測子 ABI 、変換行列の転置行列
//					ABI[nNumCluto].CorrectABI[i][j]
//						+= clust[nNumCluto].TransMatrix[0][i][k]
//						  *ABI[nNumClutplusone].PredictABI[k][l]
//						  *clust[nNumCluto].TransMatrix[0][j][l];
					// 変換行列、予測子 ABI 、変換行列の転置行列
					ABI[nNumCluto].CorrectABI[i][j]
						+= clust[nNumClutplusone].TransMatrix[0][i][k]
						  *ABI[nNumClutplusone].PredictABI[k][l]
						  *clust[nNumClutplusone].TransMatrix[0][j][l];
				}
			}
		}
	}

	// 通常の MD 及び LD の場合
/*	if (DYNMODE == MD || DYNMODE == LD)
	{*/
		for(i=0;i<6;++i)
		{
			for(j=0;j<6;++j)
			{
				// 慣性行列を足す
				ABI[nNumCluto].CorrectABI[i][j] 
				       += clust[nNumCluto].InertiaMatrix[i][j]/*kg*m^2/kg/kgm*/;
			}
		}
/*	}*/
	// BD の場合
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/*	else if (DYNMODE == BD)
	{
		for(i=0;i<6;++i)
		{
			for(j=0;j<6;++j)
			{
				// 慣性行列を足す
				ABI[nNumCluto].CorrectABI[i][j] 
				      += clust[nNumCluto].friction_tensor_tra
						*clust[nNumCluto].InertiaMatrix[i][j]
				      /*kg*m^2/kg/kgm*//*;
			}
		}
	}*/
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
}  // 0811
}

