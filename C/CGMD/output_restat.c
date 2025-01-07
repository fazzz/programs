#include <stdio.h>
#include <stdlib.h>
#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"
#include "physics.h"
#include "MD.h"

// 再スタート用のファイル出力
void output_restat(FILE *output_rest_1, 
                   FILE *output_rest_2, int nNumStep)
{
	int i;
	int alpha;
	int nNumClut;

	if ((output_rest_1 = fopen("coord.in_rest","w")) == NULL)
	{
		printf("coord.in_rest cannot open\n");
		exit(1);
	}

	if ((output_rest_2 = fopen("velo.in_rest","w")) == NULL)
	{
		printf("error velo.in_rest cannot open\n");
		exit(1);
	}

	for(i=0; i<prot.num_atom; ++i)
	{
		for(alpha=0; alpha<3; ++alpha)
		{
			fprintf(output_rest_1,"%6.3lf ",prot.coord[i][alpha]);
		}
		fprintf(output_rest_1,"\n ");
	}

 	fprintf(output_rest_1,"\n\n ");
 	fprintf(output_rest_1,"%d ",nNumStep);

	for(nNumClut=0; nNumClut<prot.DOF-1; ++nNumClut)
	{
		fprintf(output_rest_2,"%6.3lf ",clust[nNumClut].ddihedang[0]);
		fprintf(output_rest_2,"\n ");
	}

 	fprintf(output_rest_2,"\n\n ");
 	fprintf(output_rest_2,"%d ",nNumStep);
}
