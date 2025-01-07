#include <stdio.h>
#include <math.h>
#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"
#include "physics.h"
#include "MD.h"

#define nNumSpace 16

void output_file_pdb(FILE *output_c, int nNumStep)
{
	int i,ii;
	int alpha;
	int a_num;
	char *atomname_table[MAXA], nameofatom[10];

	atomname_table[1] = "C";
	atomname_table[2] = "H";
	atomname_table[3] = "O";
	atomname_table[4] = "N";
	atomname_table[5] = "S";

	if ((nNumStep%out_put_steps) == 0)
	{
		for(i=0; i<prot.num_atom; ++i)
		{
			fprintf(output_c,"ATOM    %3-d  ",i+1);

			a_num = prot.name_atom[i];
			sprintf(nameofatom,"%s",atomname_table[a_num]);

			fprintf(output_c,"%s",nameofatom);

			for(ii=0; ii<nNumSpace; ++ii)
			{
				fprintf(output_c," ");
			}

			for(alpha=0; alpha<3; ++alpha)
			{
				fprintf(output_c,"%6.3lf ",prot.coord[i][alpha]);
			}
			fprintf(output_c,"\n");
		}

		fprintf(output_c,"\n\n ");

		fclose(output_c);

	}
}
