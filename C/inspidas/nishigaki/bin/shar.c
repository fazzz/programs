/*******************************************************************************
* shar.c
* shar like file transfer from PC to UNIX
*******************************************************************************/

#include <stdio.h>
#include <string.h>

void main(argc, argv)
    int argc;
    char *argv[];
{
int i,j;
FILE *fp;
char comment[101];

        if(argc == 1){
            fprintf(stderr, "Usage : shar file1 file2 ... > transfile\n", argv[0]);
            return;
        }

        for(i=1; i < argc; i++){

            if((fp = fopen(argv[i], "r")) == NULL){
                fprintf(stderr, "Could not open %s\n", argv[i]);
		continue;
	    }

	    fprintf(stderr, "%s\n", argv[i]);
            fprintf(stdout, "sed 's/^@//' > \"%s\" << '@//E*O*F %s//'\n",
                                argv[i], argv[i]);

            while(fgets(comment, 100, fp) != NULL){
                    fputs(comment, stdout);
		    if(strrchr(comment, '\n') == NULL)
			fputs("\n", stdout);
	    }

            fprintf(stdout, "@//E*O*F %s//\n", argv[i]);

	    fclose(fp);
        }
}
