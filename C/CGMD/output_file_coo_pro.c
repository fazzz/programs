#include <stdio.h>
#include <math.h>

#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"
#include "physics.h"
#include "MD.h"

// 遂次座標データの出力ファイル
void output_file_coo_pro(FILE *output_c)
{
	int i;
	int alpha;

	for(i=0; i<prot.num_atom; ++i)
	{
		for(alpha=0; alpha<3; ++alpha)
		{
			fprintf(output_c,"%12.8lf ",prot.coord[i][alpha]);
		}
		fprintf(output_c,"\n ");
	}

 	fprintf(output_c,"\n\n ");
}
