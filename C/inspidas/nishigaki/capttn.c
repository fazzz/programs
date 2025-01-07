/*******************************************************************************
* INSPIDAS instration tool
* duplicate INCLUDE statement in source.fort files one with CMSP and the other
* CTTN header and remake INCLUDE file name in UNIX style
*******************************************************************************/

#include <stdio.h>
#include <ctype.h>

char filename[20][80] = {
	"/labo/qclib/inspidas/source/utilitydas.f",
	"That's all."};



main(argc,argv)
	int argc;
	char *argv[];
{
int i;

	for(i=0; ;i++){
		if(strcmp(filename[i],"That's all.") == 0) break;
		
		fprintf(stderr,"Working on %s\n",filename[i]);

		pass1(i);
		pass2(i);
	}
}  /* End of main() */


pass1(i)
	int i;
{
FILE *fopen();
FILE *fp_in,*fp_out;
char buff_INCLUDE[400];
int flag_INCLUDE;
char comment[121];
int j,k;

	if((fp_in = fopen(filename[i],"r")) == 0){
		printf("Error : Could not open %s\n", filename[i]);
		exit(1);
	}
	fp_out = fopen("temporaryfile","w");

	while(fgets(comment,120,fp_in) != NULL){
		/* skip Comment line */
		if(comment[0] == 'C'){
			fputs(comment,fp_out);
			continue;
		}
	
		/* search INCLUDE in comment */
		flag_INCLUDE = 0;
		for(j=6; j<7;j++){
			if(comment[j] != 'I') continue;

			if(strncmp(comment+j,"INCLUDE",7) == 0)
				flag_INCLUDE = 1;
		}

		if(!flag_INCLUDE){
			fputs(comment,fp_out);
			continue;
		}

		/* cap CTTN onto the INCLUDE statement and write the line out */
		fprintf(fp_out,"CTTN");
		fputs(comment+4,fp_out);
	}

	fclose(fp_in);
	fclose(fp_out);

}  /* End of pass1() */


pass2(i)
	int i;
{
FILE *fopen();
FILE *fp_in,*fp_out;
char buff_INCLUDE[400];
int flag_INCLUDE;
char comment[121];
int j,k;
char pathname[80],tmpstr[80];
int ptr_INCLUDE;
int flag_getpath;

	fp_in  = fopen("temporaryfile","r");
	fp_out = fopen(filename[i],"w");

	while(fgets(comment,120,fp_in) != NULL){
		fputs(comment,fp_out);
	}

}  /* End of pass2() */
