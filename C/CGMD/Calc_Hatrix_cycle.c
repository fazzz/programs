#include <stdio.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "MD.h"
#include "BD.h"

// double MatrixOfIJEle[4][4];

// 巨大行列の計算
void calc_Hatrix_cycle(void)
{
	int nNumClutI, nNumClutJ, nNumClutM;
	int alpha;

	for(nNumClutI=1;nNumClutI<=prot.DOF-1;++nNumClutI)
	{
		Hatrix[nNumClutI] = 0.0;
	}

	for(nNumClutI=1;nNumClutI<=prot.DOF-1;++nNumClutI)
	{
		for(nNumClutJ=1;nNumClutJ<=prot.DOF-1;++nNumClutJ)
		{
			for(nNumClutM=1;nNumClutM<=prot.DOF-1;++nNumClutM)
			{
				// 巨大行列の i 成分の計算
				Hatrix[nNumClutI] += sub_calc_Hatrix_cycle(nNumClutI, nNumClutJ, nNumClutM);
			}
		}
	}
}

// 巨大行列の i成分の計算
double sub_calc_Hatrix_cycle(int nNumClutI, int nNumClutJ, int nNumClutM)
{
	int alpha,i;
	int nNumClutMax;
	int nNumClut;
	int nNumClutK;

	double dummy = 0.0;
	double mat[4][4];

	nNumClutMax = nNumClutI;

	if (nNumClutMax < nNumClutJ)
	{
		nNumClutMax = nNumClutJ;
	}

	if (nNumClutMax < nNumClutM)
	{
		nNumClutMax = nNumClutM;
	}

	dummy = 0.0;
	for (nNumClutK=nNumClutMax;nNumClutK<=prot.DOF-1;++nNumClutK)
	{
		// 巨大行列の i,j 成分のひとつの項の計算
		dummy += sub_sub_calc_Hatrix_cycle(nNumClutI, nNumClutJ, nNumClutM, nNumClutK);
	}

	return dummy;
}

// 巨大行列の i,j 成分のひとつの項の計算
double sub_sub_calc_Hatrix_cycle(int nNumClutI,
							     int nNumClutJ,
							     int nNumClutM,
							     int nNumClutK)
{
	int i,j,k;
	int nNumClut;
	double matOne[4][4], matTwo[4][4], matDummy[4][4], mat[4][4];
	double Element;

	for (i=0;i<4;++i)
	{
		for (j=0;j<4;++j)
		{
			matOne[i][j] = 0.0;
			matTwo[i][j] = 0.0;
		}
	}

	calc_dot_dot_Pseduo_TransMatrix(nNumClutJ,nNumClutM,nNumClutK,matOne);

	calc_dot_Pseduo_TransMatrix(nNumClutI,nNumClutK,matTwo);

	for (i=0;i<4;++i)
	{
		for (j=0;j<4;++j)
		{
			matDummy[i][j] = 0.0;
		}
	}
	for (i=0;i<4;++i)
	{
		for (j=0;j<4;++j)
		{
			for (k=0;k<4;++k)
			{
				matDummy[i][j] += matOne[i][k]*clust[nNumClutK].PsedoInertia[k][j];
			}
		}
	}

	for (i=0;i<4;++i)
	{
		for (j=0;j<4;++j)
		{
			mat[i][j] = 0.0;
		}
	}

	for (i=0;i<4;++i)
	{
		for (j=0;j<4;++j)
		{
			for (k=0;k<4;++k)
			{
//				mat[i][j] += matDummy[i][k]*matTwo[k][j];
				mat[i][j] += matDummy[i][k]*matTwo[j][k];
			}
		}
	}

	Element = 0.0;

	for (i=0;i<4;++i)
	{
		Element += mat[i][i];
	}

	return Element*clust[nNumClutJ].ddihedang[0]*clust[nNumClutM].ddihedang[0];
}

void calc_dot_dot_Pseduo_TransMatrix(int nNumClutJ,
						             int nNumClutM,
						             int nNumClutK,
						             double mat[4][4])
{
	int i,j,k;
	int nNumClut;
	int dummy;
	double dummy_matrix[4][4], dummy_matrix2[4][4];

///////////////////////////////////////////////////////////////////////
	for (i=0;i<4;++i)
	{
		for (j=0;j<4;++j)
		{
			dummy_matrix[i][j] = ident[i][j];
		}
	}
///////////////////////////////////////////////////////////////////////

	if (nNumClutJ > nNumClutM)
	{
		dummy = nNumClutJ;
		nNumClutJ = nNumClutM;
		nNumClutM = dummy;
	}

//	if (nNumClutK >= nNumClutM && nNumClutK >= nNumClutJ /*&& nNumClutM >= nNumClutJ*/)
//	{
		for (nNumClut=1;nNumClut<=nNumClutJ;++nNumClut)
		{
			for (i=0;i<4;++i)
			{
				for (j=0;j<4;++j)
				{
					dummy_matrix2[i][j] = dummy_matrix[i][j];
				}
			}
///////////////////////////////////////////////////////////////////////
			for (i=0;i<4;++i)
			{
				for (j=0;j<4;++j)
				{
					dummy_matrix[i][j] = 0.0;
				}
			}
///////////////////////////////////////////////////////////////////////
			for (i=0;i<4;++i)
			{
				for (j=0;j<4;++j)
				{
					for (k=0;k<4;++k)
					{
						dummy_matrix[i][j] += dummy_matrix2[i][k]
											 *clust[nNumClut].PsedoTransMatrix[k][j];
					}
				}
			}
		}

		for (i=0;i<4;++i)
		{
			for (j=0;j<4;++j)
			{
				dummy_matrix2[i][j] = dummy_matrix[i][j];
			}
		}

///////////////////////////////////////////////////////////////////////
		for (i=0;i<4;++i)
		{
			for (j=0;j<4;++j)
			{
				dummy_matrix[i][j] = 0.0;
			}
		}
///////////////////////////////////////////////////////////////////////

		for (i=0;i<4;++i)
		{
			for (j=0;j<4;++j)
			{
				for (k=0;k<4;++k)
				{
					dummy_matrix[i][j] += dummy_matrix2[i][k]
										 *delta_matrix[k][j];
				}
			}
		}

		for (nNumClut=nNumClutJ+1;nNumClut<=nNumClutM;++nNumClut)
		{
			for (i=0;i<4;++i)
			{
				for (j=0;j<4;++j)
				{
					dummy_matrix2[i][j] = dummy_matrix[i][j];
				}
			}
///////////////////////////////////////////////////////////////////////
			for (i=0;i<4;++i)
			{
				for (j=0;j<4;++j)
				{
					dummy_matrix[i][j] = 0.0;
				}
			}
///////////////////////////////////////////////////////////////////////

			for (i=0;i<4;++i)
			{
				for (j=0;j<4;++j)
				{
					for (k=0;k<4;++k)
					{
						dummy_matrix[i][j] += dummy_matrix2[i][k]
											 *clust[nNumClut].PsedoTransMatrix[k][j];
					}
				}
			}
		}

		for (i=0;i<4;++i)
		{
			for (j=0;j<4;++j)
			{
				dummy_matrix2[i][j] = dummy_matrix[i][j];
			}
		}
///////////////////////////////////////////////////////////////////////
		for (i=0;i<4;++i)
		{
			for (j=0;j<4;++j)
			{
				dummy_matrix[i][j] = 0.0;
			}
		}
///////////////////////////////////////////////////////////////////////

		for (i=0;i<4;++i)
		{
			for (j=0;j<4;++j)
			{
				for (k=0;k<4;++k)
				{
					dummy_matrix[i][j] += dummy_matrix2[i][k]
										 *delta_matrix[k][j];
				}
			}
		}

		for (nNumClut=nNumClutM+1;nNumClut<=nNumClutK;++nNumClut)
		{
			for (i=0;i<4;++i)
			{
				for (j=0;j<4;++j)
				{
					dummy_matrix2[i][j] = dummy_matrix[i][j];
				}
			}
///////////////////////////////////////////////////////////////////////
			for (i=0;i<4;++i)
			{
				for (j=0;j<4;++j)
				{
					dummy_matrix[i][j] = 0.0;
				}
			}
///////////////////////////////////////////////////////////////////////

			for (i=0;i<4;++i)
			{
				for (j=0;j<4;++j)
				{
					for (k=0;k<4;++k)
					{
						dummy_matrix[i][j] += dummy_matrix2[i][k]
											 *clust[nNumClut].PsedoTransMatrix[k][j];
					}
				}
			}
		}

		for (i=0;i<4;++i)
		{
			for (j=0;j<4;++j)
			{
				mat[i][j] = dummy_matrix[i][j];
			}
		}
//	}
//	else
//	{
//		for (i=0;i<4;++i)
//		{
//			for (j=0;j<4;++j)
//			{
//				mat[i][j] = 0.0;
//			}
//		}
//	}
}

