/*******************************************************************************
* INSPIDAS instration tool
* duplicate INCLUDE statement in source.fort files one with CMSP and the other
* CTTN header and remake INCLUDE file name in UNIX style
*******************************************************************************/

#include <stdio.h>
#include <ctype.h>

char filename[20][80] = {
	"/labo/qclib/inspidas/source/dadas",
	"/labo/qclib/inspidas/source/fedrtemp",
	"/labo/qclib/inspidas/source/flxsub",
	"/labo/qclib/inspidas/source/graph",
	"/labo/qclib/inspidas/source/iolib",
	"/labo/qclib/inspidas/source/ioroutn",
	"/labo/qclib/inspidas/source/main",
	"/labo/qclib/inspidas/source/mathlib",
	"/labo/qclib/inspidas/source/mc",
	"/labo/qclib/inspidas/source/mcmain",
	"/labo/qclib/inspidas/source/minim",
	"/labo/qclib/inspidas/source/nma",
	"/labo/qclib/inspidas/source/nmanal",
	"/labo/qclib/inspidas/source/nmanal1",
	"/labo/qclib/inspidas/source/precep",
	"/labo/qclib/inspidas/source/prein",
	"/labo/qclib/inspidas/source/reg1",
	"/labo/qclib/inspidas/source/reg2v2",
	"/labo/qclib/inspidas/source/utilitie",
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

		/* write out */
		fprintf(fp_out,"CMSP");
		fputs(comment+4,fp_out);
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
		if(strncmp(comment,"CTTN",4) != 0){
			fputs(comment,fp_out);
			continue;
		}

		/* search INCLUDE in comment */
		flag_INCLUDE = 0;
		for(j=6; j<7; j++){
			if(comment[j] != 'I') continue;

			if(strncmp(comment+j,"INCLUDE",7) == 0){
				flag_INCLUDE = 1;
				break;
			}
		}

		if(!flag_INCLUDE){
			fputs(comment,fp_out);
			continue;
		}

		/* get pathname */
		ptr_INCLUDE = j;
		flag_getpath = 0;
		for(j=ptr_INCLUDE+7,k=0; j<80; j++){
			if(flag_getpath){
				if(comment[j] == ')') break;

				pathname[k] = comment[j];
				k++;
			}
			else{
				if(comment[j] == '(') flag_getpath = 1;
				continue;
			}
		}
		pathname[k] = '\0';

		/* remake pathname */
		for(k=0; ; k++){
			if(pathname[k] == '\0') break;

			pathname[k] = tolower(pathname[k]);
		}
		strcpy(tmpstr,"/labo/qclib/inspidas/common/");
		strcat(tmpstr,pathname);
		strcpy(pathname,tmpstr);

		/* write out */
/*
		strncpy(tmpstr,comment,ptr_INCLUDE);
		tmpstr[ptr_INCLUDE] = '\0';
		fprintf(fp_out,"%s",tmpstr);
		fprintf(fp_out,"INCLUDE '%s'",pathname);
		strcpy(tmpstr,comment+ptr_INCLUDE+7+j+1);
		fputs(tmpstr,fp_out);
*/
//		fprintf(fp_out,"CTTN  INCLUDE '%s'\n",pathname);
		fprintf(fp_out,"      INCLUDE '%s'\n",pathname);
	}

	fclose(fp_in);
	fclose(fp_out);

}  /* End of pass2() */
