#include <stdio.h>
#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "MD.h"
#include "BD.h"

void Jordan(int numOfRow, double Mat[100][101], double vect[100])
{
	int i,j,k,l,m,ii,jj;
	double Element[100][101];
	double ElementOfArrow;

	for (i=0;i<numOfRow;++i)
	{
		ElementOfArrow = Mat[i][i];

		for (j=0;j<numOfRow+1;++j)
		{
			Element[i][j] = Mat[i][j]/ElementOfArrow;
		}
		for (k=0;k<numOfRow;++k)
		{
			if (i != k)
			{
				for (l=0;l<numOfRow+1;++l)
				{
					Element[k][l] = Mat[k][l] - Element[i][l]*Mat[k][i];
				}
			}
		}
		for (k=0;k<numOfRow;++k)
		{
			for (l=0;l<numOfRow+1;++l)
			{
				Mat[k][l] = Element[k][l];
			}
		}

	}

	for (i=0;i<numOfRow;++i)
	{
		vect[i] = Mat[i][numOfRow];
	}

}

