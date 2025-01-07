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
FILE *fp_msp,*fp_ttn,*fp_out;
char filename[81];
char comment[81];
int i,j,k;
char msppath[81],ttnpath[81];
int ptr_slash;
static char blanks[41] = "                                        ";

	strcpy(filename,"tblmsp");
	if((fp_msp = fopen(filename,"r")) == 0){
		printf("Error : Could not open %s\n", filename);
		exit(1);
	}

	strcpy(filename,"tblttn");
	if((fp_ttn = fopen(filename,"r")) == 0){
		printf("Error : Could not open %s\n", filename);
		exit(1);
	}

	while(1){
		if(feof(fp_msp)) break;

		fgets(comment,80,fp_msp);
		sscanf(comment,"%s",msppath);
		fgets(comment,80,fp_ttn);
		sscanf(comment,"%s",ttnpath);

/*
		strcat(msppath,blanks);
		strcat(msppath,blanks);
		msppath[80] = '\0';
		strcat(ttnpath,blanks);
		strcat(ttnpath,blanks);
		ttnpath[80] = '\0';

		strcpy(comment,msppath);
		strcat(comment,ttnpath);

		fprintf(stdout,"%s\n",comment);
*/

		fprintf(stdout,"%s\n",msppath);
		fprintf(stdout,"%s\n",ttnpath);
	}

	fclose(fp_msp);
	fclose(fp_ttn);

}  /* End of main() */


