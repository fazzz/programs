#include <stdio.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"

void Initialize_Variables2(void);

void Initialize_Variables(void)
{
	int i,j;
	int nNumClut;

	for (nNumClut = 1; nNumClut < prot.DOF; ++nNumClut)
	{
		ABI[nNumClut].ABJointI = 0.0;
		zzz[nNumClut].nyu=0.0;
		zzz[nNumClut].eata=0.0;
		clust[nNumClut].dddihedang[0] = 0.0;
//		clust[nNumClut].now_deltadihedang[0] = 0.0;
		for(i=0;i<6;++i)
		{
			clust[nNumClut].sp_velo[i] = 0.0;
			clust[nNumClut].sp_acc[i] = 0.0;
			clust[nNumClut].predict_alpha[i] = 0.0;
			ABI[nNumClut].Kalman_Gain[i] = 0.0;
			zzz[nNumClut].Predictzzz[i]=0.0;
			zzz[nNumClut].Correctzzz[i]=0.0;
			clust[nNumClut].Coriolis_acc[i] = 0.0;
			clust[nNumClut].Coriolis_acc[i] = 0.0;
			/*********************************************/
                        /* clust[nNumClut].f_c.f_dihed = 0.0;	     */
                        /*********************************************/
			for(j=0;j<6;++j)
			{
				ABI[nNumClut].CorrectABI[i][j] = 0.0;
				ABI[nNumClut].PredictABI[i][j] = 0.0;
				clust[nNumClut].TransMatrix[0][i][j] = 0.0;
			}
		}
		for (i=0;i<3;++i)
		{
                        /**************************************************************/
                        /* clust[nNumClut].f_c.sp_f_clust[0].N_clust[i] = 0.0;	      */
			/* clust[nNumClut].f_c.sp_f_clust[0].f_clust[i] = 0.0;	      */
                        /**************************************************************/
			for (j=0;j<3;++j)
			{
//				clust[nNumClut].trans_A_to_CN[0][i][j] = 0.0;
			}
		}
		for(i=0;i<clust[nNumClut].num_atom_clust;++i)
		{
			for (j=0;j<3;++j)
			{
				/********************************************************/
                                /* clust[nNumClut].f_c.f_elesta[i][j] = 0.0;	        */
				/* clust[nNumClut].f_c.f_L_J[i][j] = 0.0;	        */
				/* clust[nNumClut].f_c.f_1_4_elesta[i][j] = 0.0;        */
				/* clust[nNumClut].f_c.f_1_4_L_J[i][j] = 0.0;	        */
                                /********************************************************/
			}
		}
	}
}

void Initialize_Variables2(void)
{
	int i,j;
	int nNumClut;

	for (nNumClut = 1; nNumClut < prot.DOF; ++nNumClut)
	{
		for (i=0;i<3;++i)
		{
			for (j=0;j<3;++j)
			{
				clust[nNumClut].trans_A_to_CN[0][i][j] = 0.0;
			}
		}
		for(i=0;i<clust[nNumClut].num_atom_clust;++i)
		{
			for (j=0;j<3;++j)
			{
				clust[nNumClut].xoord_clust/*[0]*/[i][j] = 0.0;
			}
		}
	}
}

