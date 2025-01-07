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
		ptr_slash = -1;
		for(i=0; ;i++){
			if(comment[i] == '\n'){
				comment[i] = ')';
				break;
			}
			else if(comment[i] == '/'){
				ptr_slash = i;
				comment[i] = '.';
			}
			else{
				comment[i] = toupper(comment[i]);
			}
		}
		if(ptr_slash != -1)
			comment[ptr_slash] = '(';

		strcpy(msppath,"INSPIDAS.");
		strcat(msppath,comment);

		fprintf(stdout,"%s\n",msppath);
	}

	fclose(fp_in);

}  /* End of main() */


