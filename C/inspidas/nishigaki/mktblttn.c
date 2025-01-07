/*******************************************************************************
* INSPIDAS instration tool
* make filename table
*******************************************************************************/

#include <stdio.h>
#include <ctype.h>



main(argc,argv)
	int argc;
	char *argv[];
{
FILE *fopen();
FILE *fp_in,*fp_out;
char filename[81];
char comment[81];
int i,j,k;
char msppath[41],ttnpath[41];
int ptr_slash;

	strcpy(filename,"tblttn");
	if((fp_in = fopen(filename,"r")) == 0){
		printf("Error : Could not open %s\n", filename);
		exit(1);
	}

	while(fgets(comment,80,fp_in) != NULL){
		strcpy(ttnpath,"/labo/qclib/inspidas/");
		strcat(ttnpath,comment);

		fputs(ttnpath,stdout);
	}

	fclose(fp_in);

}  /* End of main() */


