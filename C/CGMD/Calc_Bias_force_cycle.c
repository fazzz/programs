#include <stdio.h>
#include <math.h>
#include "f2c.h"
#include "clapack.h"

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "BD.h"
#include "MD.h"

int InvMat(double mat[6][6], double inv[6][6]);
int InvMat2(double mat[2][2], double inv[2][2]);
int InvMat3(double mat[3][3], double inv[3][3]);
int InvMat4(double mat[4][4], double inv[4][4]);

double MatHing[6];

FILE *testBF;

// Bias force ÇÃåvéZÇçsÇ§ä÷êî
void calc_Bias_force_cycle(int nNumClut){
  int nNumTerminal;
  int nNumClutTarget;
  int nNumClutOfChild[10];
  int child;
  int num_branch;

  nNumClutTarget = IndexOfABICycle[nNumClut]-1;

//	testBF = fopen("testBF.out", "a");
	// ññí[ÇÃ Bias force ÇÃåvéZ
//	if (nNumClutTarget==prot.DOF-1)
//	if (clust[nNumClutTarget].terminal == 0)
//	{
//		sub_calc_Bias_force_initial(nNumClutTarget);
//	}
//	else
//	{


  num_branch = clust[nNumClutTarget].num_branch; // 08_11
  if (nNumClutTarget==-1) num_branch=0; // 08_11
  for (child=0;child<num_branch/*clust[nNumClutTarget].num_branch**0811*/;++child){
    nNumClutOfChild[child]
      = clust[nNumClutTarget].nNumClutOfChild[child]-1;
  }
  num_branch = clust[nNumClutTarget].num_branch;
  if (nNumClutTarget==-1) num_branch=0; // 08_11
  sub_calc_Bias_force_cycle(nNumClutTarget ,
			    nNumClutOfChild,
			    num_branch);
//	}
//	fclose(testBF);

}

// ññí[ÇÃ Bias force ÇÃåvéZÇçsÇ§ä÷êî
void sub_calc_Bias_force_initial(int nNumClut){
  int i,j;
  
  int nNumClutOfChildOfThisClust[10];

  // ó\ë™ÅAèCê≥ Bias force ÇÃèâä˙âª
  //	fprintf(testBF, "Predictzzz[prot.DOF] = \n");
  for(i=0;i<6;++i){
    zzz[prot.DOF].Predictzzz[i]=0.0;
    //		fprintf(testBF, "%e \n",zzz[prot.DOF-1].Predictzzz[i]);
    //		zzz[prot.DOF-1].Correctzzz[i]=0.0;
  }
  //	fprintf(testBF,"\n");
  // eata ÇÃåvéZ
  for(i=0;i<6;++i){
    for(j=0;j<6;++j){
      clust[prot.DOF].TransMatrix[0][i][j] = 0.0;
    }
  }
  nNumClutOfChildOfThisClust[0] = prot.DOF;
  CalcCorrectzzz(nNumClut, nNumClutOfChildOfThisClust, 1);
  calceata(nNumClut);
  // nyu ÇÃåvéZ
  calcnyu(nNumClut);
  calcPredictzzz(nNumClut);
}

// Bias force ÇÃåvéZÇçsÇ§ä÷êî_1
void sub_calc_Bias_force_cycle(int nNumCluto, 
			       int nNumClutplusone[10],
			       int num_branch){
  if (nNumCluto!=-1) {
  // èCê≥éq Bias force ÇÃåvéZÇçsÇ§
  CalcCorrectzzz(nNumCluto, nNumClutplusone, num_branch);
  // eata ÇÃåvéZÇçsÇ§
  calceata(nNumCluto);
  // nyu ÇÃåvéZÇçsÇ§
  calcnyu(nNumCluto);
  // ó\ë™éq Bias force ÇÃåvéZÇçsÇ§
  calcPredictzzz(nNumCluto);
  }
}

//// Bias force ÇÃåvéZÇçsÇ§ä÷êî_2
//void sub_sub_calc_Bias_force_cycle(int nNumCluto, int nNumClutplusone)
//{
//	sub_CalcCorrectzzz(nNumCluto, nNumClutplusone);
//	calceata(nNumCluto);
//	calcnyu(nNumCluto);
//	calcPredictzzz(nNumCluto);
//}

// ó\ë™éq Bias force ÇÃåvéZÇçsÇ§ä÷êî
void calcPredictzzz(int nNumCluto){
  int i,j;
  
  // ó\ë™éq Bias force ÇÃèâä˙âªÇçsÇ§
  for(i=0;i<6;++i){
    zzz[nNumCluto].Predictzzz[i] = 0.0;
  }
  
  // KalmanGin x eata
  for(i=0;i<6;++i){
    zzz[nNumCluto].Predictzzz[i]
      = ABI[nNumCluto].Kalman_Gain[i]*zzz[nNumCluto].eata;
  }

  //	fprintf(testBF, "Predictzzz[%d] = \n", nNumCluto);
  // èCê≥éq Bias Force Çë´Ç∑
  for(i=0;i<6;++i){		
    zzz[nNumCluto].Predictzzz[i] += zzz[nNumCluto].Correctzzz[i];
    //		fprintf(testBF, "%e \n",zzz[nNumCluto].Predictzzz[i]/**1.0e20*/);
  }
//	fprintf(testBF, "\n");
//////////////////////////////////////////////////////////////
//	// ó\ë™éq ABIxcoriolis â¡ë¨ìxÇÃâ¡éZ
	for(i=0;i<6;++i)  {
	  for(j=0;j<6;++j){
//			zzz[nNumCluto].Predictzzz[i]/* N*m */
//				+=  ABI[nNumCluto].CorrectABI[i][j]
//					*clust[nNumCluto].Coriolis_acc[j]/* N*m */;
	  }
	}
//////////////////////////////////////////////////////////////
	// coriolis óÕÇÃâ¡éZ
	for(i=0;i<6;++i){
	  //			zzz[nNumCluto].Correctzzz[i]/* N*m */
	  //						 += clust[nNumCluto].Coriolis_b[i]/* N*m */;
	}
//	if (nNumCluto == prot.DOF-1)
//	{
//		for (i=0;i<6;++i)
//		{
//			zzz[nNumCluto].Predictzzz[i] = 0.0;
//		}
//	}

}

// eata ÇÃåvéZÇçsÇ§ä÷êî
void calceata(int nNumCluto){
  int i,j;
  int num_dof;
  double T;
  double vect[6];

  // eata ÇÃèâä˙âª
  zzz[nNumCluto].eata =  0.0;
  
  num_dof=/*ABI[nNumCluto].hingmat-1*/2;
  
  // ä÷êﬂãÏìÆóÕÇì«Ç›çûÇﬁ
  T /* J/rad */ = /*-(7_22)*/clust[nNumCluto].f_c.f_dihed;
  if (restflag==ON) {
    T += clust[nNumCluto].f_c.f_rest;
  }

  //  printf("%d %e\n",nNumCluto,T);

  zzz[nNumCluto].eata =  T;
  zzz[nNumCluto].eata -= zzz[nNumCluto].Correctzzz[num_dof];
  //////////////////////////////////////////////////////////////
  for(i=0;i<6;++i){
    vect[i]= 0.0;
  }

  // ó\ë™éq ABIxcoriolis â¡ë¨ìxÇÃâ¡éZ
  for(i=0;i<6;++i){
    for(j=0;j<6;++j){
      vect[i]+=  ABI[nNumCluto].CorrectABI[i][j]*clust[nNumCluto].Coriolis_acc[j];
    }
  }
//	zzz[nNumCluto].eata -=  vect[2];
//////////////////////////////////////////////////////////////
}

// n u ÇÃåvéZÇçsÇ§ä÷êî
void calcnyu(int nNumCluto){
  int num_dof,i,j,k;
  double inv[6][6],temp0[6][6],temp[3][3],temp4[4][4],inv3[3][3],inv4[4][4],temp2[2][2],inv2[2][2],add[6];
  num_dof=/*ABI[nNumCluto].hingmat-1*/2;  
  zzz[nNumCluto].nyu=zzz[nNumCluto].eata/ABI[nNumCluto].ABJointI;
  //	fprintf(testBF, "nyu[%d] = %e \n",nNumCluto , zzz[nNumCluto].nyu/**1.0e20*/);
  if(nNumCluto==0){
    for (i=0;i<6;++i)
      for (j=0;j<6;++j)
	temp0[i][j]=ABI[0].CorrectABI[i][j];
    InvMat(temp0,inv);
    for(i=0;i<6;++i) {
      acc_Term[i]=0.0;
    }
    for(i=0;i<6;++i) {
      for(j=0;j<6;++j) {
	acc_Term[i]-=inv[i][j]*zzz[0].Correctzzz[j];
      }
    }
  }

  if(nNumCluto==0 && TermMoveMode2 == 5){
    for (i=3;i<6;++i)
      for (j=3;j<6;++j)
	temp[i-3][j-3]=ABI[0].CorrectABI[i][j];
    InvMat3(temp,inv3);
    for(i=0;i<6;++i) {
      acc_Term[i]=0.0;
    }
    for(i=0;i<3;++i) {
      for(j=0;j<3;++j) {
	acc_Term[i+3]/*+*/-=inv3[i][j]*zzz[0].Correctzzz[j+3];
      }
    }
  }


  if(nNumCluto==0 && TermMoveMode2 == 6){
    for (i=2;i<6;++i)
      for (j=2;j<6;++j)
	temp4[i-2][j-2]=ABI[0].CorrectABI[i][j];
    InvMat4(temp4,inv4);
    for(i=0;i<6;++i) {
      acc_Term[i]=0.0;
    }
    for(i=0;i<4;++i) {
      for(j=0;j<4;++j) {
	acc_Term[i+2]/*+*/-=inv4[i][j]*zzz[0].Correctzzz[j+2];
      }
    }
  }

  if(nNumCluto==0 && TermMoveMode2 == 7){
    for (i=0;i<3;++i)
      for (j=0;j<3;++j)
	temp[i][j]=ABI[0].CorrectABI[i][j];
    InvMat3(temp,inv3);
    for(i=0;i<6;++i) {
      acc_Term[i]=0.0;
    }
    for(i=0;i<3;++i) {
      for(j=0;j<3;++j) {
	acc_Term[i]/*+*/-=inv3[i][j]*zzz[0].Correctzzz[j];
      }
    }
  }

  if(nNumCluto==0 && TermMoveMode2 == 8){
    zzz[0].nyu=-zzz[0].Correctzzz[2]/ABI[0].CorrectABI[2][2];
    for(i=0;i<6;++i) {
      acc_Term[i]=0.0;
    }
    acc_Term[2]=zzz[0].nyu;
  }

  if(nNumCluto==0 && TermMoveMode2 == 9){
    for (i=2;i<4;++i)
      for (j=2;j<4;++j)
	temp2[i-2][j-2]=ABI[0].CorrectABI[i][j];
    InvMat2(temp2,inv2);
    for(i=0;i<6;++i) {
      acc_Term[i]=0.0;
    }
    for(i=0;i<2;++i) {
      for(j=0;j<2;++j) {
	acc_Term[i+2]/*+*/-=inv2[i][j]*zzz[0].Correctzzz[j+2];
      }
    }
  }

  if(nNumCluto==0 && TermMoveMode2 == 10){
    for (i=0;i<6;++i)
      for (j=0;j<6;++j)
	temp0[i][j]=ABI[0].CorrectABI[i][j];
    InvMat(temp0,inv);
    for(i=0;i<6;++i) {
      acc_Term[i]=0.0;
      acc_Term2[i]=0.0;
    }
    add[0]=clust[0].InertiaMatrix[5][5]*
      (  (-vel_Term[2]*clust[0].qCOM[2]-vel_Term[1]*clust[0].qCOM[1])*vel_Term[3]
        +vel_Term[1]*clust[0].qCOM[0]*vel_Term[4]
        +vel_Term[2]*clust[0].qCOM[0]*vel_Term[5]);
    add[1]=clust[0].InertiaMatrix[5][5]*
      (  vel_Term[0]*clust[0].qCOM[1]*vel_Term[3]
        +(-vel_Term[2]*clust[0].qCOM[2]-vel_Term[0]*clust[0].qCOM[0])*vel_Term[4]
        +vel_Term[2]*clust[0].qCOM[1]*vel_Term[5]);
    add[2]=clust[0].InertiaMatrix[5][5]*
      (  vel_Term[0]*clust[0].qCOM[2]*vel_Term[3]
        +vel_Term[1]*clust[0].qCOM[2]*vel_Term[4]
        +(-vel_Term[0]*clust[0].qCOM[0]-vel_Term[1]*clust[0].qCOM[1])*vel_Term[5]);

    add[3]=clust[0].InertiaMatrix[5][5]*(-vel_Term[2]*vel_Term[3+1]+vel_Term[1]*vel_Term[3+2]);
    add[4]=clust[0].InertiaMatrix[5][5]*( vel_Term[2]*vel_Term[3+0]-vel_Term[0]*vel_Term[3+2]);
    add[5]=clust[0].InertiaMatrix[5][5]*(-vel_Term[1]*vel_Term[3+0]+vel_Term[0]*vel_Term[3+1]);

    for(i=0;i<6;++i) {
      for(j=0;j<6;++j) {
	acc_Term[i]-=inv[i][j]*(zzz[0].Correctzzz[j]+add[j]);
      }
    }
    for(i=0;i<6;++i) {
      for(j=0;j<6;++j) {
	acc_Term2[i]-=inv[i][j]*(zzz[0].Correctzzz[j]);
      }
    }
  }

  if(nNumCluto==0 && TermMoveMode2==11){
    for (i=0;i<6;++i)
      for (j=0;j<6;++j)
	temp0[i][j]=ABI[0].CorrectABI[i][j];
    InvMat(temp0,inv);
    for(i=0;i<6;++i) {
      acc_Term[i]=0.0;
    }
    add[0]=0.0;
    add[1]=0.0;
    add[2]=0.0;
    add[3]=(-vel_Term[2]*vel_Term[3+1]+vel_Term[1]*vel_Term[3+2]);
    add[4]=( vel_Term[2]*vel_Term[3+0]-vel_Term[0]*vel_Term[3+2]);
    add[5]=(-vel_Term[1]*vel_Term[3+0]+vel_Term[0]*vel_Term[3+1]);

    for(i=0;i<6;++i) {
      for(j=0;j<6;++j) {
	acc_Term[i]-=inv[i][j]*zzz[0].Correctzzz[j];
      }
      acc_Term[i]-=add[i];
    }
  }

  if(nNumCluto==0 && TermMoveMode2==12){
    for (i=0;i<6;++i)
      for (j=0;j<6;++j)
	temp0[i][j]=ABI[0].CorrectABI[i][j];
    InvMat(temp0,inv);
    for(i=0;i<6;++i) {
      acc_Term[i]=0.0;
      acc_Term2[i]=0.0;
      acc_Term3[i]=0.0;
    }
    add[0]=0.0;
    add[1]=0.0;
    add[2]=0.0;
    add[3]=(-vel_Term[2]*vel_Term[3+1]+vel_Term[1]*vel_Term[3+2]);
    add[4]=( vel_Term[2]*vel_Term[3+0]-vel_Term[0]*vel_Term[3+2]);
    add[5]=(-vel_Term[1]*vel_Term[3+0]+vel_Term[0]*vel_Term[3+1]);

    for(i=0;i<6;++i) {
      for(j=0;j<6;++j) {
	acc_Term2[i]-=inv[i][j]*zzz[0].Correctzzz[j];
      }
    }
    for(i=0;i<6;++i) {
      acc_Term[i]=acc_Term2[i]-add[i];
    }
    for (i=0;i<3;++i) {
      for (j=0;j<3;++j) {
	acc_Term3[i]+=Rot_Term[j][i]*acc_Term2[j];
      }
      for (j=0;j<3;++j) {
	acc_Term3[i+3]+=Rot_Term[i][j]*acc_Term2[j+3];
	//	acc_Term3[i+3]+=Rot_Term[i][j]*acc_Term[j];
      }
    }
  }

}

// èCê≥éq Bias force ÇÃåvéZÇçsÇ§ä÷êî_1
void CalcCorrectzzz(int nNumCluto, 
		    int nNumClutplusone[10], 
		    int num_branch){
  int i,j;
  int child, nNumNow;

  // èCê≥éq Bias force ÇÃèâä˙âª
  for(i=0;i<6;++i) {
    zzz[nNumCluto].Correctzzz[i] = 0.0;
  }

  for (child=0;child<num_branch;++child){
    nNumNow = nNumClutplusone[child];
    if (nNumNow != -2){
///////////////////////////////////////////////////////////////////////////////
      // ïœä∑çsóÒÇÃèÊéZ
      for(i=0;i<6;++i){
	for(j=0;j<6;++j){
	  zzz[nNumCluto].Correctzzz[i]
	    +=   clust[nNumNow].TransMatrix[0][i][j]
	    * zzz[nNumNow].Predictzzz[j];
	}
      }
    }
///////////////////////////////////////////////////////////////////////////////
  }

////	if (DYNMODE == MD || DYNMODE == LD)
////	if (DYNMODE == MD || DYNMODE == LD)
////	{
//			// ó\ë™éq ABIxcoriolis â¡ë¨ìxÇÃâ¡éZ
  for(i=0;i<6;++i) {
    for(j=0;j<6;++j){
      ////				zzz[nNumCluto].Correctzzz[i]/* N*m */
      ////						 +=   ABI[nNumCluto].PredictABI[i][j]
      ////							 *clust[nNumCluto].Coriolis_acc[j]/* N*m */;
      zzz[nNumCluto].Correctzzz[i]/* N*m */
	+=   ABI[nNumCluto].CorrectABI[i][j]
	*clust[nNumCluto].Coriolis_acc[j]/* N*m */;
    }
  }
////	}

//		// coriolis óÕÇÃâ¡éZ
  for(i=0;i<6;++i){
    zzz[nNumCluto].Correctzzz[i]/* N*m */
      += clust[nNumCluto].Coriolis_b[i]/* N*m */;
  }

  // óÕÇÃå∏éZ_1
  for(i=0;i<3;++i){
    zzz[nNumCluto].Correctzzz[i] /* N */
      /*+*/-= clust[nNumCluto].f_c.sp_f_clust[0].N_clust[i]/* N */;
  }

//	fprintf(testBF, "Correctzzz[%d] = \n", nNumCluto);
  // óÕÇÃå∏éZ_2
  for(i=3;i<6;++i){
    zzz[nNumCluto].Correctzzz[i] /* N*m */
      /*+*/-= clust[nNumCluto].f_c.sp_f_clust[0].f_clust[i-3]/* N*m */;
  }
//	for(i=0;i<6;++i)
//	{
//		fprintf(testBF, "%e \n",zzz[nNumCluto].Correctzzz[i]/**1.0e20*/);
//	}

}

// èCê≥éq Bias force ÇÃåvéZÇçsÇ§ä÷êî_2
//void sub_CalcCorrectzzz(int nNumCluto, int nNumClutplusone)
//{
//	int i,j;
//	int alpha, alpha2;
//
//	double TransMatrix[6][6];
//	double Tensorxvelo[6];

//	// èCê≥éq Bias force ÇÃèâä˙âª
//	for(i=0;i<6;++i)
//	{
//		zzz[nNumCluto].Correctzzz[i] = 0.0;
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
//			TransMatrix[i][j] = -1*clust[nNumCluto].TransMatrix[0][i][j];
//		}
//	}
//
//	for (i=0;i<6;++i)
//	{
//		TransMatrix[i][i] += 2;
//	}
//
//	// ïœä∑çsóÒÇÃèÊéZ
//	for(i=0;i<6;++i)
//	{
//		for(j=0;j<6;++j)
//		{
//			zzz[nNumCluto].Correctzzz[i]
//					 +=   TransMatrix[i][j]
//					    * zzz[nNumClutplusone].Predictzzz[j];
//		}
//	}
//
//	// ó\ë™éq ABIxcoriolis â¡ë¨ìxÇÃâ¡éZ
//	for(i=0;i<6;++i)
//	{
//		for(j=0;j<6;++j)
//		{
//			zzz[nNumCluto].Correctzzz[i]/* N*m */
//					 +=   ABI[nNumCluto].PredictABI[i][j]
//						 *clust[nNumCluto].Coriolis_acc[j]/* N*m */;
//		}
//	}
//
//	// coriolis óÕÇÃâ¡éZ
//	for(i=0;i<6;++i)
//	{
//		zzz[nNumCluto].Correctzzz[i]/* N*m */
//					 += clust[nNumCluto].Coriolis_b[i]/* N*m */;
//	}
//
//	// óÕÇÃå∏éZ_1
//	for(i=0;i<3;++i)
//	{
////		zzz[nNumCluto].Correctzzz[i] /* N */
////		             -= clust[nNumCluto].f_c.sp_f_clust[0].N_clust[i]/* N */;
//	}
//
//	// óÕÇÃå∏éZ_2
//	for(i=3;i<6;++i)
//	{
////		zzz[nNumCluto].Correctzzz[i] /* N*m */
////		         -= clust[nNumCluto].f_c.sp_f_clust[0].f_clust[i-3]/* N*m */;
//	}
//
//	if (DYNMODE == LD)
//	{
//
//		// ñÄéCóÕÇÃâ¡éZ
//		for(i=0;i<6;++i)
//		{
//			for(j=0;j<3;++j)
//			{
//				zzz[nNumCluto].Correctzzz[i]
//			             += clust[nNumCluto].InertiaMatrix[i][j]
//						   	*clust[nNumCluto].friction_tensor_tra
//						    *clust[nNumCluto].sp_velo[j];
//			}
//		}
//
//// BD
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//		// óhìÆóÕÇÃå∏éZ
//		for(i=0;i<6;++i)
//		{
//			zzz[nNumCluto].Correctzzz[i]
//		         	-= clust[nNumCluto].f_c.f_Brownian[i];
//		}
//	}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//}

int InvMat(double mat[6][6], double inv[6][6]) {

  int i,j,k;
  double A[36];
  double test[6][6];
  static long int m=6,n=6,lda=6,info,piv[6],lwork=6;
  static double work[6];

  k=0;
  for(i=0;i<6;++i) {
    for(j=0;j<6;++j) {
      A[k]=mat[j][i];
      ++k;
    }
  }

  dgetrf_(&m,&n,A,&lda,piv,&info);

  dgetri_(&n,A,&lda,piv,work,&lwork,&info);

  k=0;
  for(i=0;i<6;++i) {
    for(j=0;j<6;++j) {
      inv[j][i]=A[k];
      ++k;
    }
  }

  for(i=0;i<6;++i) {
    for(j=0;j<6;++j) {
      test[i][j] = 0.0;
    }
  }

  for(i=0;i<6;++i) {
    for(j=0;j<6;++j) {
      for(k=0;k<6;++k) {
	test[i][j]+=mat[i][k]*inv[k][j];
      }
    }
  }


  return 0;

}

int InvMat3(double mat[3][3], double inv[3][3]) {

  int i,j,k;
  double A[9];
  double test[3][3];
  static long int m=3,n=3,lda=3,info,piv[3],lwork=3;
  static double work[3];

  k=0;
  for(i=0;i<3;++i) {
    for(j=0;j<3;++j) {
      A[k]=mat[j][i];
      ++k;
    }
  }

  dgetrf_(&m,&n,A,&lda,piv,&info);

  dgetri_(&n,A,&lda,piv,work,&lwork,&info);

  k=0;
  for(i=0;i<3;++i) {
    for(j=0;j<3;++j) {
      inv[j][i]=A[k];
      ++k;
    }
  }

  for(i=0;i<3;++i) {
    for(j=0;j<3;++j) {
      test[i][j] = 0.0;
    }
  }

  for(i=0;i<3;++i) {
    for(j=0;j<3;++j) {
      for(k=0;k<3;++k) {
	test[i][j]+=mat[i][k]*inv[k][j];
      }
    }
  }


  return 0;

}

int InvMat4(double mat[4][4], double inv[4][4]) {

  int i,j,k;
  double A[16];
  double test[4][4];
  static long int m=4,n=4,lda=4,info,piv[4],lwork=4;
  static double work[4];

  k=0;
  for(i=0;i<4;++i) {
    for(j=0;j<4;++j) {
      A[k]=mat[j][i];
      ++k;
    }
  }

  dgetrf_(&m,&n,A,&lda,piv,&info);

  dgetri_(&n,A,&lda,piv,work,&lwork,&info);

  k=0;
  for(i=0;i<4;++i) {
    for(j=0;j<4;++j) {
      inv[j][i]=A[k];
      ++k;
    }
  }

  for(i=0;i<4;++i) {
    for(j=0;j<4;++j) {
      test[i][j] = 0.0;
    }
  }

  for(i=0;i<4;++i) {
    for(j=0;j<4;++j) {
      for(k=0;k<4;++k) {
	test[i][j]+=mat[i][k]*inv[k][j];
      }
    }
  }


  return 0;

}

int InvMat2(double mat[2][2], double inv[2][2]) {

  int i,j,k;
  double A[4];
  double test[2][2];
  static long int m=2,n=2,lda=2,info,piv[2],lwork=2;
  static double work[2];

  k=0;
  for(i=0;i<2;++i) {
    for(j=0;j<2;++j) {
      A[k]=mat[j][i];
      ++k;
    }
  }

  dgetrf_(&m,&n,A,&lda,piv,&info);

  dgetri_(&n,A,&lda,piv,work,&lwork,&info);

  k=0;
  for(i=0;i<2;++i) {
    for(j=0;j<2;++j) {
      inv[j][i]=A[k];
      ++k;
    }
  }

  for(i=0;i<2;++i) {
    for(j=0;j<2;++j) {
      test[i][j] = 0.0;
    }
  }

  for(i=0;i<2;++i) {
    for(j=0;j<2;++j) {
      for(k=0;k<2;++k) {
	test[i][j]+=mat[i][k]*inv[k][j];
      }
    }
  }


  return 0;

}
