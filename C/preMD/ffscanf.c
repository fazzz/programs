#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

//#include "Vis_MD.h"

int ffscanf(FILE *filename, /*const char *stream,*/ char type);

int ffscanf(FILE *filename, /*const char *stream,*/ char type)
{
	int i;
	char line;

	 fscanf(filename, "%s", &type);
	if (*type == '#')
	{
		while ((i=fgetc(filename))!='\0' || (i=fgetc(filename))!='\n')
		{
			;
		}
	}
	else
	{
		return (int)type;
	}
	return (int)type;
}

int main()
{
	FILE *input;
	char l;
	if ((input = fopen("stdres.dat","r")) == NULL)
	{
		printf("stdres.dat ;cannot open\n");
		exit(1);
	}
	ffscanf(input, &l);
	fclose(input);

	printf("%s", l);
	return 0;
}
