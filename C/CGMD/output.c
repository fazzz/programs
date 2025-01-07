#include <stdio.h>
#include <math.h>

#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"
#include "physics.h"
#include "MD.h"

// ���X�e�b�v�̍\�����o��
void output(int nNumStep)
{
	int nNumAtom;
	int alpha;

	// ����̃X�e�b�v���Ƃɏo��
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
