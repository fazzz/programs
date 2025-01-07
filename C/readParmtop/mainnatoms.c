#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "PT.h"

void usage(char *progname);

#define ON 1
#define OFF 0

int main(int argc, char *argv[]) {
  int i,j;
  int outflag=OFF,nflag=OFF,resflag=OFF;
  int natom,ires;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *parmfilename,*outputfilename,*progname;
  FILE *parmfile,*outputfile;

  while((c=getopt(argc,argv,"hnrf:"))!=-1) {
    switch(c) {
    case 'f':
      outflag=ON;
      outputfilename=optarg;
      break;
    case 'n':
      nflag=ON;
      break;
    case 'r':
      resflag=ON;
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 1) {
    USAGE(progname);
    exit(1);
  }
  parmfilename = *argv;

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);

  natom=AP.NATOM;

  if (outflag==OFF) {
    ires=0;
    for (i=0;i<natom;++i) {
      if ( nflag==OFF) 
	printf("%s ",AP.IGRAPH[i]);
      else 
	for (j=0;j<4;++j) 
	  printf("%c",AP.IGRAPH[i][j]);
      if (AP.IPRES[ires+1]<i)
	ires+=1;
      if (resflag==ON)
	if ( nflag==OFF) 
	  printf("%s ",AP.LABERES[ires]);
	else 
	  for (j=0;j<4;++j) 
	    printf("%c",AP.LABERES[ires][j]);
      if (nflag==ON)
	printf("\n");
    }
    printf("\n");
  }
  /************************************************/
  /* else {					  */
  /*   outputfile=efopen(outputfilename,"w");	  */
  /*   for (i=0;i<natom;++i) 			  */
  /*     fprintf(outputfile,"%s ",AP.IGRAPH[i]);  */
  /*   if (resflag==ON)				  */
  /*     fprintf(outputfile,"%s ",AP.LABERES[i]); */
  /*   fprintf(outputfile,"\n");		  */
  /*   fclose(outputfile);			  */
  /* }						  */
  /************************************************/
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-r] res name \n");
  printf("[-f outputfilename  ] outputfilemode \n");
  printf("[-h] help \n");
  printf("%s [-f outputfilename ]  [-h] parmfilename \n",progname);
}
