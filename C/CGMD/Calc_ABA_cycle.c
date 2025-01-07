#include <stdio.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "MD.h"
#include "BD.h"

double MatHing[6];
double MatHingTranspose[6];

//FILE *outtest;

// Articulated-Body Inertia �A
// Kalman Gain �̌v�Z���s���֐�
void calc_ABA_cycle(int nNumClut)
{
	int nNumTerminal;
	int nNumClutTarget;
	int nNumClutOfChild[10];
	int child;
	nNumClutTarget = IndexOfABICycle[nNumClut]- 1;

//	outtest = fopen("testACI.out", "a");
	// ���[�� ABI �AKG �̌v�Z
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

// ���[�� ABI �AKG �̌v�Z���s���֐�
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

// ABI �̌v�Z���s���֐�
void sub_calc_ABI_cycle(int nNumCluto, // ���݂̍��̂̃C���f�b�N�X
	                    int nNumClutplusone[10], // ���݂̍��̂̐e���̂̃C���f�b�N�X
                        int num_branch)
{
	// �C���q ABI �̌v�Z���s��
	CalcCorrectABI(nNumCluto, nNumClutplusone, num_branch);
	// Joint Inertia �̌v�Z���s��
	CalcD(nNumCluto);
	// KG �̌v�Z���s��
	calcKalmanGain(nNumCluto);
	// �\���q ABI �̌v�Z���s��
	calcPredictABI(nNumCluto);
}

// ABI �̌v�Z���s���֐�
void sub_sub_calc_ABI_cycle(int nNumCluto, int nNumClutplusone)
{
	// �C���q ABI �̌v�Z���s��
	sub_CalcCorrectABI(nNumCluto, nNumClutplusone);
	// Joint Inertia �̌v�Z���s��
	CalcD(nNumCluto);
	// KG �̌v�Z���s��
	calcKalmanGain(nNumCluto);
	// �\���q ABI �̌v�Z���s��
	calcPredictABI(nNumCluto);
}

// �\���q ABI �̌v�Z���s���֐�///////////////////////////////////////
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
	// �\���q ABI �̏�����
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
	// BD �̏ꍇ
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

//	// tau �̌v�Z
//	for (i=0;i<6;++i)
//	{
//		for (j=0;j<6;++j)
//		{
//			tau[i][j] =   ident[i][j]
//                        - ABI[nNumCluto].Kalman_Gain[i]
//                         *MatHing[j];
//		}
//	}

	// �\���q ABI �̌v�Z
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

// KG �̌v�Z���s���֐�
void calcKalmanGain(int nNumCluto){
  int i,j;
  int num_dof;
  // KG �̏�����
  if (nNumCluto!=-1) { // 0811
  for(i=0;i<6;++i){
    ABI[nNumCluto].Kalman_Gain[i] = 0.0;
  }
  num_dof=/*ABI[nNumCluto].hingmat-1*/2;
  // KG �̌v�Z���s��
  for(i=0; i<6; ++i) {
    ABI[nNumCluto].Kalman_Gain[i]= ABI[nNumCluto].CorrectABI[i][num_dof]/ABI[nNumCluto].ABJointI;
  }
  }
}

// Joint Inertia �̌v�Z���s���֐�
void CalcD(int nNumCluto){
  int num_dof;
  if (nNumCluto!=-1) { // 0811
  num_dof=/*ABI[nNumCluto].hingmat-1*/2;
  ABI[nNumCluto].ABJointI = ABI[nNumCluto].CorrectABI[num_dof][num_dof];
  }
}

// �C���q ABI �̌v�Z���s���֐�_1
void CalcCorrectABI(int nNumCluto, 
		    int nNumClutplusone[10], 
		    int num_branch){
  int i,j,k,l;
  int child,nNumNow;
  double mat[6][6];

  // �C���q ABI �̏�����
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
	      // �ϊ��s��A�\���q ABI �A�ϊ��s��̓]�u�s��
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
//					// �ϊ��s��A�\���q ABI �A�ϊ��s��̓]�u�s��
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
//					// �ϊ��s��A�\���q ABI �A�ϊ��s��̓]�u�s��
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
//					// �ϊ��s��A�\���q ABI �A�ϊ��s��̓]�u�s��
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
      // �����s��𑫂�
      ABI[nNumCluto].CorrectABI[i][j] 
	+= clust[nNumCluto].InertiaMatrix[i][j]
	/*kg*m^2/kg/kgm*/;
      //	fprintf(outtest,"%e ",ABI[nNumCluto].CorrectABI[i][j]/**1.0e40*/);
    }
    //	fprintf(outtest,"\n");
  }
  }
//	fprintf(outtest,"\n");
	// BD �̏ꍇ
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/*	else if(DYNMODE == BD)
	{
		for(i=0;i<6;++i)
		{
			for(j=0;j<6;++j)
			{
				// �����s��𑫂�
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

// �C���q ABI �̌v�Z���s���֐�_2
void sub_CalcCorrectABI(int nNumCluto, int nNumClutplusone)
{
	int i,j,k,l;

	double TransMatrix[6][6];

	if (nNumCluto != -1) { // 0811 

	// �C���q Bias force �̏�����
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

	// �C���q ABI �̏�����
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
//					// �ϊ��s��A�\���q ABI �A�ϊ��s��̓]�u�s��
//					ABI[nNumCluto].CorrectABI[i][j]
//						+= clust[nNumCluto].TransMatrix[0][i][k]
//						  *ABI[nNumClutplusone].PredictABI[k][l]
//						  *clust[nNumCluto].TransMatrix[0][j][l];
					// �ϊ��s��A�\���q ABI �A�ϊ��s��̓]�u�s��
					ABI[nNumCluto].CorrectABI[i][j]
						+= clust[nNumClutplusone].TransMatrix[0][i][k]
						  *ABI[nNumClutplusone].PredictABI[k][l]
						  *clust[nNumClutplusone].TransMatrix[0][j][l];
				}
			}
		}
	}

	// �ʏ�� MD �y�� LD �̏ꍇ
/*	if (DYNMODE == MD || DYNMODE == LD)
	{*/
		for(i=0;i<6;++i)
		{
			for(j=0;j<6;++j)
			{
				// �����s��𑫂�
				ABI[nNumCluto].CorrectABI[i][j] 
				       += clust[nNumCluto].InertiaMatrix[i][j]/*kg*m^2/kg/kgm*/;
			}
		}
/*	}*/
	// BD �̏ꍇ
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/*	else if (DYNMODE == BD)
	{
		for(i=0;i<6;++i)
		{
			for(j=0;j<6;++j)
			{
				// �����s��𑫂�
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

