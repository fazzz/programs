#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string.h>
#include <ctype.h>

#include "EF.h"

#define ON 1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,ii,jj,k;
  int flag1=OFF,flag2=OFF,aflag=OFF,wrappedflag=OFF;
  int c1,c2;
  int numlen1,numlen2;
  char num1[20],num2[20];
  double a,b,d;

  char *progname;
  char *inpfile1name,*inpfile2name,*outfilename;
  FILE *inpfile1,*inpfile2,*outfile;

  char *line1,*line2,*dummy;
  char ch;
  size_t len1=0,len2=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"hag"))!=-1) {
    switch(c) {
    case 'a':
      aflag=ON;
      break;
    case 'g':
      wrappedflag=ON;
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  argc-=optind;
  argv+=optind;

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  inpfile1name  = *argv;
  inpfile2name = *++argv;
  outfilename = *++argv;

  inpfile1=efopen(inpfile1name,"r");
  inpfile2=efopen(inpfile2name,"r");
  outfile=efopen(outfilename,"w");

  while ( (c1=getline(&line1,&len1,inpfile1))!=-1 && (c2=getline(&line2,&len2,inpfile2))!=-1 ) {
    i=0;
    ii=0;
    while (isspace(line1[i]))
      ++i;
    while (isspace(line2[ii]))
      ++ii;
    numlen1=strlen(line1);
    numlen2=strlen(line2);
    while ( i < numlen1 || ii < numlen2  ) {
      while (isspace(line1[i]) && i < numlen1) {
	++i;
      }
      while (isspace(line2[ii]) && ii < numlen2 ) {
	++ii;
      }
      flag1=OFF;
      j=0;
      for (k=0;k<20;++k)
	num1[k]='\0';
      while (!isspace(line1[i])  && i < numlen1 ) {
	num1[j]=line1[i];
	++i;
	++j;
	flag1=ON;
      }
      a=atof(num1);
      jj=0;
      flag2=OFF;
      for (k=0;k<20;++k)
	num2[k]='\0';
      while (!isspace(line2[ii]) && ii < numlen2 ) {
	num2[jj]=line2[ii];
	++ii;
	++jj;
	flag2=ON;
      }
      if ( flag1 == ON && flag2 == ON ) {
	b=atof(num2);
	d=a-b;
	if ( aflag == ON ) {
	  //	  d=abs(d);
	  if ( d < 0.0   )
	    d=-1.0*d;
	}
	if ( wrappedflag == ON ) {
	  //	  d=abs(d);
	  if ( d < 0.0   )
	    d=-1.0*d;
	  if ( d > 180.0   )
	    d=180.0*2.0-d;
	  if ( d < 0.0)
	    d=-1.0*d;
	}
	fprintf(outfile,"%20.16e ",d);
      }
    }
    fprintf(outfile,"\n");
  }


  fclose(inpfile1);
  fclose(inpfile2);
  fclose(outfile);

  return 0;
}

void USAGE(char *progname) {
  printf("-a -- absolute value\n");
  printf("-g -- wrapped value [-180:180]\n");
  printf("-h -- help\n");
  printf("USAGE: %s inpfile2name inpfile2name outfilename \n", progname);
}

