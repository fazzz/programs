/*******************************************************************************
* INSPIDAS instration tool
* dasarc.c
* make archived INSPIDAS file
*******************************************************************************/

#include <stdio.h>
#include <fcntl.h>

struct{
	char dir[81];
	char file[40][81];
	int num_file;
} msppath[40];
int num_msppath;
char ttnpath[40][81];



main(argc,argv)
	int argc;
	char *argv[];
{
FILE *fopen();
FILE *fp_in,*fp_out;
char filename[81];
char comment[121],comment2[121];
int i,j,k,l;
char dir[81],dir_prev[81],file[81];
int flag_getdir;
int ptr_msppath,num_file;

	/* read filename table and fill msppath struct */
	strcpy(filename,"tbl");
	if((fp_in = fopen(filename,"r")) == 0){
		printf("Error : Could not open %s\n", filename);
		exit(1);
	}

	dir_prev[0] = '\0';
	for(i=0; ;i++){
		if(feof(fp_in)) break;

		fgets(comment,80,fp_in);
		fgets(comment2,80,fp_in);
		sscanf(comment2,"%s",ttnpath[i]);

		/* get directory name and filename on MSP */
		flag_getdir = 1;
		for(j=0,k=0,l=0; ;j++){
			if(flag_getdir){
				if(comment[j] == '('){
					dir[k] = '\0';
					flag_getdir = 0;
					continue;
				}
				else{
					dir[k++] = comment[j];
				}
			}
			else{
				if(comment[j] == ')'){
					file[l] = '\0';
					break;
				}
				else{
					file[l++] = comment[j];
				}
			}
		}

		if(strcmp(dir,dir_prev) != 0){
			num_msppath++;
			strcpy(dir_prev,dir);
		}

		ptr_msppath = num_msppath-1;
		strcpy(msppath[ptr_msppath].dir,dir);
		num_file = msppath[ptr_msppath].num_file++;
		strcpy(msppath[ptr_msppath].file[num_file],file);
	}

	fclose(fp_in);

	/* write out job stream */
	fp_out = fopen("PACKDAS","w");

	fprintf(fp_out,"//PACKDAS JOB E50214,#,CLASS=A,NOTIFY=E50214\n");

	for(i=0,k=0; i<num_msppath; i++){
		strcpy(dir,msppath[i].dir);
		num_file = msppath[i].num_file;
		fprintf(fp_out,"/* %s\n",dir);
		fprintf(fp_out,"//  EXEC UPDATE,PARM=NEW\n");
		fprintf(fp_out,"//SYSIN DD *\n");

		for(j=0; j<num_file; j++){
			strcpy(file,msppath[i].file[j]);
			fprintf(fp_out,"./ ADD NAME=%s\n",file);
			fprintf(fp_out,"./ NUMBER NEW1=10,INCR=10\n");

			strcpy(filename,ttnpath[k++]);
			if((fp_in = fopen(filename,"r")) == 0){
				fprintf(stderr,"Error : Could not open %s\n", filename);
				continue;
			}
			while(fgets(comment,120,fp_in) != NULL){
				fputs(comment,fp_out);
			}
			fclose(fp_in);
		}

		fprintf(fp_out,"//SYSUT2 DD DSN=%s,UNIT=PUB,\n",dir);
		fprintf(fp_out,"//  SPACE=(TRK,(100,10,40),RLSE),\n");
		fprintf(fp_out,"//  DISP=(NEW,CATLG),\n");
		fprintf(fp_out,"//  DCB=(RECFM=FB,LRECL=80,BLKSIZE=6320)\n");
	}

	fprintf(fp_out,"//\n");

	fclose(fp_out);

}  /* End of main() */

