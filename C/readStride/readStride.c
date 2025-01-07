
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "readStride.h"

#include "PDB.h"
#include "EF.h"
#include "PT.h"

#define ON 1
#define OFF 0

int readStride(FILE *stridefile,SD *SDdata,int numres) {
  int i,j,c,n,flag;
  double f,pmflag;
  char type[3];

  char *line,dummy;
  size_t len=0;

  i=-1;
  n=-1;

  while ((c=getc(stridefile))!=-1){
    if (c=='\n')
      n=-1;
    else
      ++n;
    if (n >= 0 && n < 3 ) type[n]=c;
    if (n==4) {
      flag=OFF;
      if (strncmp(type,"ASG",3)==0) {
	++i;
	flag=ON;
      }
    }

    j=0;
    if (n >= 5 && n < 8 ) {
      SDdata[i].resname[j]=c;
      ++j;
    }

    if (n == 24 ) {
      SDdata[i].OL_type_2nd_st=c;
    }

    j=0;
    if (n >= 30 && n < 39 ) {
      SDdata[i].type_2nd_st[j]=c;
      ++j;
    }

    pmflag=1.0;
    if (39<= n && n <49) {
      if (c==' ')
	f=0;
      else if (c=='-')
	pmflag=-1.0;
      else if (isdigit(c))
	f=(double)(c-'0')*pow(10.0,(49-n-1));
      SDdata[i].phi+=f;
    }
    SDdata[i].phi=pmflag*SDdata[i].phi;

    if (31<= n && n <39) {
      SDdata[i].psi+=f;
    }
    SDdata[i].psi=pmflag*SDdata[i].psi;

    pmflag=1.0;
    if (49<= n && n <59) {
      if (c==' ')
	f=0;
      else if (c=='-')
	pmflag=-1.0;
      else if (isdigit(c))
	f=(double)(c-'0')*pow(10.0,(49-n-1));
      SDdata[i].type_2nd_st[j]+=c;
    }
    SDdata[i].psi=pmflag*SDdata[i].psi;

    pmflag=1.0;
    if (59<= n && n <69) {
      if (c==' ')
	f=0;
      else if (c=='-')
	pmflag=-1.0;
      else if (isdigit(c))
	f=(double)(c-'0')*pow(10.0,(49-n-1));
      SDdata[i].area+=f;
    }
    SDdata[i].area=pmflag*SDdata[i].area;
  }

  fclose(stridefile);
}
