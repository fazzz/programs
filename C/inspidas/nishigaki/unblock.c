
/*******************************************************************************
* INSPIDAS instration tool
* unblock.c
* dd conv=unblock
*******************************************************************************/

#include <stdio.h>
#include <string.h>



main(argc,argv)
	int argc;
	char *argv[];
{
FILE *fopen();
FILE *fp_in,*fp_out;
char filename[81];
char comment[121];
int i,j,k;

	if(argc == 1){
		fprintf(stderr,"Usage : unblock filename\n");
		exit(0);
	}
	
	strcpy(filename,argv[1]);
	if((fp_in = fopen(filename,"r")) == 0){
		printf("Error : Could not open %s\n", filename);
		exit(1);
	}
	fp_out = fopen("unblockedfile","w");

	while(fgets(comment,120,fp_in) != NULL){
		for(i=strlen(comment)-2; i>=0; i--){
			if(comment[i] == ' ') continue;
			
			comment[i+1] = '\n';
			comment[i+2] = '\0';

			break;
		}
		fputs(comment,fp_out);
	}
}  /* End of main() */
