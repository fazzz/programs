#include <stdio.h>
#include <math.h>

#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"
#include "physics.h"
#include "MD.h"

// 現ステップの構造を出力
void output(int nNumStep)
{
	int nNumAtom;
	int alpha;

	// 所定のステップごとに出力
	if ((nNumStep % out_put_steps) == 0)
	{
		printf("prot.coord = \n");

		for(nNumAtom=0; nNumAtom<prot.num_atom; ++nNumAtom)
		{
			for(alpha=0; alpha<3; ++alpha)
			{
				printf("%6.3lf ",prot.coord[nNumAtom][alpha]);
			}
			printf("\n");
		}
	}
}
