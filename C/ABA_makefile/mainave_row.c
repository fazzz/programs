#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string.h>
#include <ctype.h>

#include "EF.h"

#define ON 1
#define OFF 0

#define INTEGER 0
#define DECIMAL 1
#define WRITE   2
#define SPACE   3
#define ADDLINE 4
#define ADDLINEWRITE 5
#define EXPO    6

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,argnum;
  int sign,n,ne;
  double f,v,d;
  int e;
  int signe=1.0;
  int flag,vflag=ON,oflag=OFF,hflag=OFF,iflag=OFF;
  
  int status;

  int nl,nc,numlen,ncmax,nf;
  char num[20];
  double *ave,*var,*index;
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
  while((c=getopt(argc,argv,"hvoei"))!=-1) {
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
    case 'i':
      iflag=ON;
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
  index=(double *)gcemalloc(sizeof(double)*1);

  inpfile=efopen(inpfilename,"r");
  if (hflag==ON) getline(&line,&len,inpfile);

  nl=0;
  ncmax=1;
  status=SPACE;
  d=0.0;
  f=0.0;
  nf=0;
  e=0;
  n=0;
  sign=1.0;
  while ( (c=getc(inpfile))!=-1 ) {
    i=0;
    nc=0;
    if(isspace(c) && c!='\n') {
      if (status==INTEGER || status==DECIMAL || status==EXPO ) {
	status=WRITE;
      }
      else {
	status=SPACE;
      }
    }
    else if (isdigit(c)) {
      if (status==INTEGER)
	d=d*10+c-'0';
      else if (status==DECIMAL) {	
	++nf;
	f=f+(c-'0')*pow(0.1,nf);
      }
      else if (status==EXPO) {	
	e=e*10+c-'0';
      }
      else if (status==SPACE || status==WRITE) {
	status=INTEGER;
	d=c-'0';
      }	
    }
    else if (c=='.') {
      if (status!=EXPO) status=DECIMAL;
      else printf("error\n");
    }
    else if (c=='+') {
      if (status==INTEGER || status==DECIMAL ) {
	status=WRITE;
      }
      else if (status==EXPO) {
	signe=1;
      }
      else {
	status=INTEGER;
	sign=1.0;
      }
    }
    else if (c=='-') {
      if (status==INTEGER || status==DECIMAL ) {
	status=WRITE;
      }
      else if (status==EXPO) {
	signe=-1;
      }
      else {
	status=INTEGER;
	sign=-1.0;
      }
    }
    else if (c=='e' || c=='E') {
      status=EXPO;
    }
    else if (c=='\n') {
      if (status==INTEGER || status==DECIMAL || status==EXPO ) {
	status=ADDLINEWRITE;
      }
      else 
	status=ADDLINE;
    }
    if (status==WRITE || status==ADDLINEWRITE) {
      e=signe*e;
      v=sign*(f+d);
      v=v*pow(10.0,e);
      d=0.0;
      f=0.0;
      nf=0.0;
      e=0;
      signe=1.0;
      sign=1.0;
      ++n;
      if (iflag==OFF) {
	ave[nl]=((n-1)*ave[nl]+v)/n;
	if ( vflag == ON ) var[nl]=((n-1)*var[nl]+v*v)/n;
      }
      else {
	if (n>1) {
	  ave[nl]=((n-2)*ave[nl]+v)/(n-1);
	  if ( vflag == ON ) var[nl]=((n-2)*var[nl]+v*v)/(n-1);
	}
	else index[nl]=v;
      }
    }
    if (status==ADDLINEWRITE || status==ADDLINE) {
      n=0;
      ++nl;
      status=SPACE;
      ave=(double *)gcerealloc(ave,sizeof(double)*nl);
      if ( vflag == ON ) var=(double *)gcerealloc(var,sizeof(double)*nl);
      if ( iflag == ON ) index=(double *)gcerealloc(index,sizeof(double)*nl);
    }
  }

  if (oflag==ON) {
    if (iflag==OFF) {
      for (i=0;i<nl;++i) {
	printf("%20.16e ",ave[i]);
	if ( vflag == ON )
	  printf("%20.16e ",sqrt(var[i]-ave[i]*ave[i]));
	printf("\n");
      }
    }
    else  {
      for (i=0;i<nl;++i) {
	printf("%20.16e ",index[i]);
	printf("%20.16e ",ave[i]);
	if ( vflag == ON )
	  printf("%20.16e ",sqrt(var[i]-ave[i]*ave[i]));
	printf("\n");
      }
    }
  }
  else {
    outfile1=efopen(outfilename1,"w");
    if ( vflag == ON ) outfile2=efopen(outfilename2,"w");
    for (i=0;i<nl;++i) fprintf(outfile1,"%20.16e ",ave[i]);
    if ( vflag == ON )
      for (i=0;i<nl;++i) fprintf(outfile2,"%20.16e ",sqrt(var[i]-ave[i]*ave[i]));
    fclose(inpfile);
    fclose(outfile1);
    if ( vflag == ON ) fclose(outfile2);
  }

  return 0;
}

void USAGE(char *progname) {
  printf("-o -- ouputtype\n");
  printf("-v -- wo variance calc.\n");
  printf("-h -- help\n");
  printf("USAGE: %s inpfilename outfilename(ave) outfilename2(var) \n", progname);
}

