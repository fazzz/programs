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
#define EXP 4

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,argnum;
  int sign,n;
  double f,v,d,fe,signe,fse;
  int flag,vflag=ON,oflag=OFF,hflag=OFF;
  
  int status;

  int nl,nc,numlen,ncmax,nf;
  char num[20];
  double *ave,*var;
  double a;

  char *progname;
  char *inpfilename,*outfilename1,*outfilename2;
  FILE *inpfile,*outfile1,*outfile2;

  char *line;
  char ch;
  size_t len=0;

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
  status=SPACE;
  d=0.0;
  f=0.0;
  nf=0;
  n=0;
  sign=1.0;
  signe=1.0;
  while ( (c=getc(inpfile))!=-1 ) {
    i=0;
    nc=0;
    if(isspace(c) && c!='\n') {
      if (status==INTEGER || status==DECIMAL || status==EXP ) {
	++n;
	v=sign*(f+d);
	fe=signe*fe;
	fse=v*pow(10,fe);
	d=0.0;
	f=0.0;
	fe=0.0;
	nf=0;
	sign=1.0;
	signe=1.0;
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
      else if (status==EXP) {	
	fe=fe*10+c-'0';
      }
      else if (status==SPACE || status==WRITE) {
	status=INTEGER;
	d=c-'0';
      }	
    }
    else if (c=='.') {
      status=DECIMAL;
    }
    else if (c=='e') {
      status=EXP;
      fe=0.0;
      signe=1.0;
    }
    else if (c=='+') {
      if (status==INTEGER || status==DECIMAL ) {
	++n;
	v=sign*(f+d);
	d=0.0;
	f=0.0;
	nf=0.0;
	fe=0.0;
	sign=1.0;
	signe=1.0;
	status=WRITE;
      }
      else if (status!=EXP) {
	status=INTEGER;
      }

      if (status!=EXP) {
	sign=1.0;
      }
      else {
	signe=1.0;
      }
    }
    else if (c=='-') {
      if (status==INTEGER || status==DECIMAL ) {
	++n;
	v=sign*(f+d);
	d=0.0;
	f=0.0;
	nf=0.0;
	fe=0.0;
	sign=1.0;
	signe=1.0;
	status=WRITE;
      }
      else if (status!=EXP) {
	status=INTEGER;
      }

      if (status!=EXP) {
	sign=-1.0;
      }
      else {
	signe=-1.0;
      }
    }
    else if (c=='\n') {
      if (n>ncmax) ncmax=n;
      n=0;
      ++nl;
      status=SPACE;
    }
    if (status==WRITE) {
      if (nl==0) {
	ave=(double *)gcerealloc(ave,sizeof(double)*n);
	if ( vflag == ON )
	  var=(double *)gcerealloc(var,sizeof(double)*n);
      }
      ave[n-1]=(nl*ave[n-1]+fse/*v*/)/(nl+1);
      //      ave[n-1]=(nl*ave[n-1]+fse/*v*/)/(nl+1);
      if ( vflag == ON )
	var[n-1]=(nl*var[n-1]+fse*fse/*v*v*/)/(nl+1);
	//	var[n-1]=(nl*var[n-1]+(fse-ave[n-1])*(fse-ave[n-1]))/(nl+1);
    }
  }

  if (oflag==ON) {
    for (i=0;i<ncmax;++i) {
      printf("%20.16e ",ave[i]);
      if ( vflag == ON )
      	printf("%20.16e ",sqrt(var[i]-ave[i]*ave[i]));
      //	printf("%20.16e ",sqrt(var[i]));
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
  printf("-h -- help\n");
  printf("USAGE: %s inpfilename outfilename(ave) outfilename2(var) \n", progname);
}


