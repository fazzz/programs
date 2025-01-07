/*******************************************************************************
* fcut.c
* cut out the specified region from a very large file
*******************************************************************************/

#include <stdio.h>
#include <string.h>


main(argc,argv)
	int argc;
	char *argv[];
{
FILE *fopen();
FILE *fp;
int i,line_from,line_to;
char filen[80],comment[1025];	

	if(argc == 1){
		fprintf(stderr, "Usage : fcut line_from line_to filename\n");
		exit(0);
	}

/*
	fscanf(stdin, "%d %d %s", &line_from,&line_to,filen);
*/

	line_from = atoi(argv[1]);
	line_to   = atoi(argv[2]);
	strcpy(filen, argv[3]);

	if(!(fp = fopen(filen, "r"))){
		printf("Error : Could not open %s\n", filen);
		exit(0);
	}
	
	/* skip */
	for(i=1; i<line_from; i++)
		fgets(comment, 1024, fp);
		
	/* cut out and write */
	for(i=line_from; i<=line_to; i++){
		fgets(comment, 1024, fp);
		fputs(comment, stdout);
	}

	fclose(fp);
} /* End of main() */

