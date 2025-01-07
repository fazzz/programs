
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "readAOUT.h"

#define ON  1
#define OFF 0

#define FMODE 0
#define EMODE 1
#define GMODE 2

double *readOUT(FILE *inputfile,char *type,int *numstep) {
  int i,j,k;
  int c;
  int flag1=OFF,flag2;
  double *data;
  double f,g,s,e;

  data=(double *)gcemalloc(sizeof(double)*1);

  i=0;
  j=0;
  while ((c=getc(inputfile)) != EOF) {
    if (flag1==OFF) {
      if (c==type[i]) {
	++i;
      }
      else {
	i=0;
      }
      if (type[i]=='\0') {
	i=0;
	k=1;
	flag1=ON;
	flag2=FMODE;
	f=0.0;
	s=1.0;
	g=0.0;
	e=0.0;
      }
    }
    else {
      if (c >= '0' && c <= '9') {
	if (flag2==FMODE)
	  f=f*10.0+c-'0';
	else if (flag2==GMODE) {
	  g=g+(c-'0')/pow(10,k);
	  ++k;
	}
	else if (flag2==EMODE) {
	  e=e*10.0+c-'0';
	}
      }
      else if (c=='.') {
	flag2=GMODE;
      }
      else if (c == '-') {
	s=-1.0;
      }
      else if (c == 'e' || c == 'E') {
	getc(inputfile);
	flag2=EMODE;
      }
      else if (c == '\b' || c == '\t' || c == ' ' || c == '\n') {
	;
      }
      else {
	flag1=OFF;
      }

      if (flag1==OFF) {
	f=s*(f+g)*pow(10.0,e);
	data=(double *)gcerealloc(data,sizeof(double)*(j+1));
	data[j]=f;
	++j;
      }      
    }
  }

  *numstep=j;

  return data;
}

double *readOUTwnum(FILE *inputfile,char *type,int numstep) {
  int i,j,k;
  int c;
  int flag1=OFF,flag2;
  double *data;
  double f,g,s,e;

  data=(double *)gcemalloc(sizeof(double)*1);

  i=0;
  j=0;
  while ((c=getc(inputfile)) != EOF) {
    if (flag1==OFF) {
      if (c==type[i]) {
	++i;
      }
      else {
	i=0;
      }
      if (type[i]=='\0') {
	i=0;
	k=1;
	flag1=ON;
	flag2=FMODE;
	f=0.0;
	s=1.0;
	g=0.0;
	e=0.0;
      }
    }
    else {
      if (c >= '0' && c <= '9') {
	if (flag2==FMODE)
	  f=f*10.0+c-'0';
	else if (flag2==GMODE) {
	  g=g+(c-'0')/pow(10,k);
	  ++k;
	}
	else if (flag2==EMODE) {
	  e=e*10.0+c-'0';
	}
      }
      else if (c=='.') {
	flag2=GMODE;
      }
      else if (c == '-') {
	s=-1.0;
      }
      else if (c == 'e' || c == 'E') {
	getc(inputfile);
	flag2=EMODE;
      }
      else if (c == '\b' || c == '\t' || c == ' ' || c == '\n') {
	;
      }
      else {
	flag1=OFF;
      }

      if (flag1==OFF) {
	f=s*(f+g)*pow(10.0,e);
	data=(double *)gcerealloc(data,sizeof(double)*(j+1));
	data[j]=f;
	++j;
	if (j>numstep)
	  c=EOF;
      }      
    }
  }

  if (j<numstep) {
    printf("error about inputfile file\n");
    exit(1);
  }

  return data;
}

double *readOUTrow(FILE *inputfile,char *type,int numstep) {
  int i,j,k;
  int c;
  int flag1=OFF,flag2;
  double *data;
  double f,g,s,e;

  data=(double *)gcemalloc(sizeof(double)*1);

}

double *readOUTwab(FILE *inputfile,char *type,int *numstep,int numabondon) {
  int i,j,k;
  int nl=0;
  int c;
  int flag1=OFF,flag2;
  double *data;
  double f,g,s,e;

  data=(double *)gcemalloc(sizeof(double)*1);

  i=0;
  j=0;
  while ((c=getc(inputfile)) != EOF) {
    if (flag1==OFF) {
      if (c==type[i]) {
	++i;
      }
      else {
	i=0;
      }
      if (type[i]=='\0') {
	i=0;
	k=1;
	flag1=ON;
	flag2=FMODE;
	f=0.0;
	s=1.0;
	g=0.0;
	e=0.0;
      }
    }
    else {
      if (c >= '0' && c <= '9') {
	if (flag2==FMODE)
	  f=f*10.0+c-'0';
	else if (flag2==GMODE) {
	  g=g+(c-'0')/pow(10,k);
	  ++k;
	}
	else if (flag2==EMODE) {
	  e=e*10.0+c-'0';
	}
      }
      else if (c=='.') {
	flag2=GMODE;
      }
      else if (c == '-') {
	s=-1.0;
      }
      else if (c == 'e' || c == 'E') {
	getc(inputfile);
	flag2=EMODE;
      }
      else if (c == '\b' || c == '\t' || c == ' ' || c == '\n') {
	;
      }
      else {
	flag1=OFF;
      }

      if (flag1==OFF) {
	if (nl>numabondon) {
	  f=s*(f+g)*pow(10.0,e);
	  data=(double *)gcerealloc(data,sizeof(double)*(j+1));
	  data[j]=f;
	  ++j;
	}
	++nl;
      }      
    }
  }

  *numstep=j;

  return data;
}

double *readOUTwabwnum(FILE *inputfile,char *type,int numstep,int numabondon) {
  int i,j,k;
  int nl=0;
  int c;
  int flag1=OFF,flag2;
  double *data;
  double f,g,s,e;

  data=(double *)gcemalloc(sizeof(double)*1);

  i=0;
  j=0;
  while ((c=getc(inputfile)) != EOF) {
    if (flag1==OFF) {
      if (c==type[i]) {
	++i;
      }
      else {
	i=0;
      }
      if (type[i]=='\0') {
	i=0;
	k=1;
	flag1=ON;
	flag2=FMODE;
	f=0.0;
	s=1.0;
	g=0.0;
	e=0.0;
      }
    }
    else {
      if (c >= '0' && c <= '9') {
	if (flag2==FMODE)
	  f=f*10.0+c-'0';
	else if (flag2==GMODE) {
	  g=g+(c-'0')/pow(10,k);
	  ++k;
	}
	else if (flag2==EMODE) {
	  e=e*10.0+c-'0';
	}
      }
      else if (c=='.') {
	flag2=GMODE;
      }
      else if (c == '-') {
	s=-1.0;
      }
      else if (c == 'e' || c == 'E') {
	getc(inputfile);
	flag2=EMODE;
      }
      else if (c == '\b' || c == '\t' || c == ' ' || c == '\n') {
	;
      }
      else {
	flag1=OFF;
      }

      if (flag1==OFF) {
	if ( nl > numabondon) {
	  f=s*(f+g)*pow(10.0,e);
	  data=(double *)gcerealloc(data,sizeof(double)*(j+1));
	  data[j]=f;
	  ++j;
	}      
	if (nl>numstep)
	  c=EOF;
	++nl;
      }
    }
  }

  if (j<numstep-numabondon-1) {
    printf("error about inputfile file\n");
    exit(1);
  }

  return data;
}
