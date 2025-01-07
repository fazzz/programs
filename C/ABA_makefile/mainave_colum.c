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
  int i,j,argnum;
  int flag,vflag=ON,oflag=OFF,hflag=OFF;
  int nl,nc,numlen,ncmax;
  char num[20];
  double *ave,*var;
  double a;

  char *progname;
  char *inpfilename,*outfilename1,*outfilename2;
  FILE *inpfile,*outfile1,*outfile2;

  char *line;
  char ch;
  size_t len;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"hvoe"))!=-1) {
    switch(c) {
    case 'v':
      vflag=OFF;
      break;
    case 'o':
      oflag=ON;
      break;
    case 'e':
      hflag=ON;
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

  if (oflag==ON) argnum=1;
  if (oflag==OFF) argnum=2;
  if (argc < argnum) {
    USAGE(progname);
    exit(1);
  }
  inpfilename  = *argv;
  outfilename1 = *++argv;
  if ( vflag == ON )
    outfilename2 = *++argv;

  ave=(double *)gcemalloc(sizeof(double)*1);
  var=(double *)gcemalloc(sizeof(double)*1);

  inpfile=efopen(inpfilename,"r");
  if (hflag==ON) getline(&line,&len,inpfile);

  nl=0;
  ncmax=1;
  while ( (c=getline(&line,&len,inpfile))!=-1 ) {
    i=0;
    nc=0;
    c=isspace(line[i]);
    while (c!=0)
      ++i;
    numlen=strlen(line);
    while ( i < numlen ) {
      while (isspace(line[i]) && i < numlen) {
	++i;
      }
      flag=OFF;
      j=0;
      while (!isspace(line[i])  && i < numlen ) {
	num[j]=line[i];
	++i;
	++j;
	flag=ON;
      }
      a=atof(num);
      if ( flag == ON  ) {
	++nc;
	if (nl==0) {
	  ave=(double *)gcerealloc(ave,sizeof(double)*nc);
	  if ( vflag == ON )
	    var=(double *)gcerealloc(var,sizeof(double)*nc);
	}
	if (nc > ncmax)
	  ncmax=nc;
	ave[nc-1]=(nl*ave[nc-1]+a)/(nl+1);
	if ( vflag == ON )
	  var[nc-1]=(nl*var[nc-1]+a*a)/(nl+1);
      }
    }
    ++nl;
  }

  if (oflag==ON) {
    for (i=0;i<ncmax;++i) {
      printf("%20.16e ",ave[i]);
      if ( vflag == ON )
	printf("%20.16e ",sqrt(var[i]-ave[i]*ave[i]));
      printf("\n");
    }
  }
  else {
    outfile1=efopen(outfilename1,"w");
    if ( vflag == ON ) outfile2=efopen(outfilename2,"w");
    for (i=0;i<ncmax;++i) fprintf(outfile1,"%20.16e ",ave[i]);
    if ( vflag == ON )
      for (i=0;i<ncmax;++i) fprintf(outfile2,"%20.16e ",sqrt(var[i]-ave[i]*ave[i]));
    fclose(inpfile);
    fclose(outfile1);
    if ( vflag == ON ) fclose(outfile2);
  }

  return 0;
}

void USAGE(char *progname) {
  printf("-o -- ouputtype\n");
  printf("-v -- wo variance calc.\n");
  printf("-e -- headerflag\n");
  printf("-h -- help\n");
  printf("USAGE: %s inpfilename outfilename(ave) outfilename2(var) \n", progname);
}

