#include <stdio.h>
#include <math.h>
#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"
#include "physics.h"

void output_file_initial_pdb(void)
{
	int i,ii,j;
	int a_num;
	char *atomname_table[MAXA], nameofatom[10];

	FILE *output_c;

	atomname_table[1] = "C";
	atomname_table[2] = "H";
	atomname_table[3] = "O";
	atomname_table[4] = "N";
	atomname_table[5] = "S";
	
	if ((output_c = fopen("coo_pro_initial.pdb", "w")) == NULL)
	{
		printf("error coo_pro_initial.pdb cannot open\n");
		exit(1);
	}

	for(i=0; i<prot.num_atom; ++i)
	{
		fprintf(output_c,"ATOM    %3-d  ",i);

		a_num = prot.name_atom[i];
		sprintf(nameofatom,"%s",atomname_table[a_num]);

		fprintf(output_c,"%s",nameofatom);

		for(ii=0; ii<16; ++ii)
		{
			fprintf(output_c," ");
		}

		for(j=0; j<3; ++j)
		{
			fprintf(output_c,"%6.3lf ",prot.coord[i][j]);
		}
		fprintf(output_c,"\n");
	}

	fprintf(output_c,"\n\n ");

	fclose(output_c);
}
