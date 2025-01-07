
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "OPT.h"

int getopt(int flag[],int numflag,int argc,char *argv[],char *inputfilename[], char *USAGE) {
  int i,j;
  int flag2=1;
  char *option;
  int opt;

  if (argc > 1){
    i = 1;
    while (++i < argc){
      if ( (*++argv)[0] == '-'){
	option = *argv;
	opt = *++option;
	flag2=1;
	for (j=0;j<numflag;++j) {
	  if (opt == flag[j]){
	    ++i;
	    if (*++argv != NULL){
	      inputfilename[j] = *argv;
	      flag2=0;
	    }
	    else {
	      flag2=1;
	    }
	  }
	}
	if (flag2==1) {
	  printf("%s\n",USAGE);
	  exit(1);
	}
      }
    }
  }
}

  
