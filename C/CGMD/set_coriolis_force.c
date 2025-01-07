#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "physics.h"
#include "MD.h"

// Coriolis Force の設定の補助を行う関数
void sub_set_coriolis_force(int nNumClut,
	                        int nNumClutminousone);


// Coriolis Force の設定を行う関数
void set_coriolis_force(int nNumClut)
{

	int i;
	int nNumClutOfParent;

	nNumClutOfParent = clust[nNumClut].nNumClutOfParent-1;

//	// クラスタが末端のとき
//	if (nNumClutOfParent == 0)
//	{
//		for(i=0;i<6;++i)
//		{
//			clust[nNumClut].Coriolis_b[i] = 0.0;
//		}
//	}
	// クラスタが末端のとき
	if (nNumClut == 0) {
          if (TermMoveMode==OFF) {
	    for(i=0;i<6;++i) {
	      clust[nNumClut].Coriolis_b[i] = 0.0;
	    }
          }
          else {
	    sub_set_coriolis_force(nNumClut,-1);
	  }
	}
	else {
	  sub_set_coriolis_force(nNumClut, nNumClutOfParent);
	}

}

// Coriolis Force の設定の補助を行う関数
void sub_set_coriolis_force(int nNumClut,
	                        int nNumClutminousone)
{
  int i,j,k;
	int alpha, alpha2,alpha3;

	double omega[3];
	double Inertia[3][3];

	double omegaproduct[3][3];
	double MasproductS[3];
	double matOmegaOmega[3][3];
	double InertiaProductOmega[3];
	double MasproductSProductOmega[3];

	double omegaqCOMvelo[3];
	double omegaqCOM[3][3];
	double qCOMproduct[3][3];


	for (i=0;i<6;++i)
	{
		clust[nNumClut].Coriolis_b[i] = 0.0;
	}

	for (alpha=0;alpha<3;++alpha)
	{
		InertiaProductOmega[alpha] = 0.0;
		MasproductSProductOmega[alpha] = 0.0;
	}

//上半分の設定///////////////////////////////////////////////////////
//	omegaproduct[0][0]= 0.0;
//	omegaproduct[1][0]=-clust[nNumClut].sp_velo[2]/*m*/;
//	omegaproduct[2][0]= clust[nNumClut].sp_velo[1]/*m*/;
//	omegaproduct[0][1]= clust[nNumClut].sp_velo[2]/*m*/;
//	omegaproduct[1][1]= 0.0;
//	omegaproduct[2][1]=-clust[nNumClut].sp_velo[0]/*m*/;
//	omegaproduct[0][2]=-clust[nNumClut].sp_velo[1]/*m*/;
//	omegaproduct[1][2]= clust[nNumClut].sp_velo[0]/*m*/;
//	omegaproduct[2][2]= 0.0;

	omegaproduct[0][0]= 0.0;
	omegaproduct[0][1]=-clust[nNumClut].sp_velo[2]/*m*/;
	omegaproduct[0][2]= clust[nNumClut].sp_velo[1]/*m*/;
	omegaproduct[1][0]= clust[nNumClut].sp_velo[2]/*m*/;
	omegaproduct[1][1]= 0.0;
	omegaproduct[1][2]=-clust[nNumClut].sp_velo[0]/*m*/;
	omegaproduct[2][0]=-clust[nNumClut].sp_velo[1]/*m*/;
	omegaproduct[2][1]= clust[nNumClut].sp_velo[0]/*m*/;
	omegaproduct[2][2]= 0.0;
//	omegaproduct[0][0]= 0.0;
//	omegaproduct[0][1]=-clust[nNumClut].ddihedang[0]/*m*/;
//	omegaproduct[0][2]= 0.0;
//	omegaproduct[1][0]= clust[nNumClut].ddihedang[0]/*m*/;
//	omegaproduct[1][1]= 0.0;
//	omegaproduct[1][2]= 0.0;
//	omegaproduct[2][0]= 0.0;
//	omegaproduct[2][1]= 0.0;
//	omegaproduct[2][2]= 0.0;

	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
//			InertiaProductOmega[alpha]
//			            += 2.0*clust[nNumClut].Inertia_clust[alpha][alpha2]
//			              *clust[nNumClut].sp_velo[alpha2];
			InertiaProductOmega[alpha]
			            += clust[nNumClut].InertiaMatrix[alpha][alpha2]
			              *clust[nNumClut].sp_velo[alpha2];
		}
	}

	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
//			clust[nNumClut].Coriolis_b[alpha]
//			            += omegaproduct[alpha][alpha2]
//			              *InertiaProductOmega[alpha2];
////			clust[nNumClut].Coriolis_b[alpha]
////			            -= omegaproduct[alpha][alpha2]
////			              *InertiaProductOmega[alpha2];
			clust[nNumClut].Coriolis_b[alpha]
			            +=/*2.0**/omegaproduct[alpha][alpha2]
			              *InertiaProductOmega[alpha2];
		}
	}
/////////////////////////////////////////////////////////////////////

//下半分の設定///////////////////////////////////////////////////////

//	MasproductS[0] = clust[nNumClut].InertiaMatrix[4][2];
//	MasproductS[1] = clust[nNumClut].InertiaMatrix[5][0];
//	MasproductS[2] = clust[nNumClut].InertiaMatrix[3][1];
	MasproductS[0] = clust[nNumClut].qCOM[0]*Sum_Mass(nNumClut);
	MasproductS[1] = clust[nNumClut].qCOM[1]*Sum_Mass(nNumClut);
	MasproductS[2] = clust[nNumClut].qCOM[2]*Sum_Mass(nNumClut);
//	MasproductS[0] = -clust[nNumClut].qCOM[0]*Sum_Mass(nNumClut);
//	MasproductS[1] = -clust[nNumClut].qCOM[1]*Sum_Mass(nNumClut);
//	MasproductS[2] = -clust[nNumClut].qCOM[2]*Sum_Mass(nNumClut);

	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			MasproductSProductOmega[alpha]
			            += omegaproduct[alpha][alpha2]
			              *MasproductS[alpha2];
		}
	}

	for (alpha=0;alpha<3;++alpha)\
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
//			clust[nNumClut].Coriolis_b[alpha+3]
//			            += /*0.5**/omegaproduct[alpha][alpha2]
//			              *MasproductSProductOmega[alpha2];
//			clust[nNumClut].Coriolis_b[alpha+3]
//			            -= omegaproduct[alpha][alpha2]
//			              *MasproductSProductOmega[alpha2];
		}
	}
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			matOmegaOmega[alpha][alpha2] = 0.0;
		}
	}

	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			for (alpha3=0;alpha3<3;++alpha3)
			{
				matOmegaOmega[alpha][alpha2] += omegaproduct[alpha][alpha3]
											   *omegaproduct[alpha3][alpha2];
			}
		}
	}
	for (alpha=0;alpha<3;++alpha)
	{
		for (alpha2=0;alpha2<3;++alpha2)
		{
			clust[nNumClut].Coriolis_b[alpha+3]
			             += matOmegaOmega[alpha][alpha2]
			              *MasproductS[alpha2];
		}
	}
	if (TermMoveMode2 == 4) {
	  //	  if (nNumClut==0) {
	  qCOMproduct[0][0]= 0.0;
	  qCOMproduct[0][1]=-clust[nNumClut].qCOM[2]*Sum_Mass(nNumClut);
	  qCOMproduct[0][2]= clust[nNumClut].qCOM[1]*Sum_Mass(nNumClut);
	  qCOMproduct[1][0]= clust[nNumClut].qCOM[2]*Sum_Mass(nNumClut);
	  qCOMproduct[1][1]= 0.0;
	  qCOMproduct[1][2]=-clust[nNumClut].qCOM[0]*Sum_Mass(nNumClut);
	  qCOMproduct[2][0]=-clust[nNumClut].qCOM[1]*Sum_Mass(nNumClut);
	  qCOMproduct[2][1]= clust[nNumClut].qCOM[0]*Sum_Mass(nNumClut);
	  qCOMproduct[2][2]= 0.0;

	  for (i=0;i<3;++i)
	    for (j=0;j<3;++j)
	      omegaqCOM[i][j]=0.0;
	  for (i=0;i<3;++i)
	    for (j=0;j<3;++j)
	      for (k=0;k<3;++k)
		omegaqCOM[i][j]+=omegaproduct[i][k]*qCOMproduct[k][j];

	  for (i=0;i<3;++i)
	      omegaqCOMvelo[i]=0.0;
	  for (i=0;i<3;++i)
	    for (j=0;j<3;++j)
	      omegaqCOMvelo[i]+=omegaqCOM[i][j]*clust[nNumClut].sp_velo[j];

          for (i=0;i<3;++i)
	    clust[nNumClut].Coriolis_b[i]+=omegaqCOMvelo[i];

	  for (alpha=0;alpha<3;++alpha)
	    clust[nNumClut].Coriolis_b[alpha+3]=0.0;
	  //	  }
	}
/////////////////////////////////////////////////////////////////////


}
