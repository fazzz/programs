#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "PDB.h"
#include "EF.h"

int readPDB(char *filename) {  
  int i,j,c,n,na,d,ATOMDATAFLAG=0;
  double f,pmflag;
  FILE *pdbfile;

  pdbfile=efopen(filename,"r");

  n=-1;
  na=0;
  while ((c=getc(pdbfile))!=-1/*EOF*/){
    if (c=='\n') {
      n=-1;
      if (ATOMDATAFLAG==1) {
	++na;
	PD.serial[na]=0;
	PD.resSeq[na]=0;
	for (i=0;i<3;++i) {
	  PD.coord[na][i]=0.0;
	}
	ATOMDATAFLAG==0;
      }
    }
    else {
      ++n;
    }
    if (n < 6) {
      PD.recordname[na][n]=c;
    }
    PD.recordname[na][4]=0;
    if (n==6) {
      if (strcmp(PD.recordname[na],"ATOM")==0){
	ATOMDATAFLAG=1;
      }
    }
    if (ATOMDATAFLAG==1){
      if (6 <= n && n <11) {
	if (c==' ') {
	  d=0;
	}
	else if (isdigit(c)) {
	  d=(c-'0')*pow(10,(11-n-1));
	}
	PD.serial[na]+=d;
      }
      if (11<= n && n <16)
	PD.name[na][n-11]=c;
      if (n == 16)
	PD.altLOC[na]=c;
      if (17<= n && n <20)
	PD.resname[na][n-17]=c;
      if (n == 21)
	PD.ChainID[na]=c;
      if (22<= n && n <26) {
	if (c==' ') {
	  d=0;
	}
	else if (isdigit(c)) {
	  d=(c-'0')*pow(10,(26-n-1));
	}
	PD.resSeq[na]+=d;
      }
      if (n==26)
	PD.iCode[na]=c;
      if (n==30)
	pmflag=1.0;
      if (30<= n && n <34) {
	if (c==' ') {
	  f=0;
	}
	else if (c=='-') {
	  pmflag=-1.0;
	}
	else if (isdigit(c)) {
	  f=(double)(c-'0')*pow(10.0,(34-n-1));
	}
	PD.coord[na][0]+=f;
      }
      if (35<= n && n <38) {
	if (c==' ') {
	  f=0;
	}
	else if (isdigit(c)) {
	  f=(double)(c-'0')*pow(10.0,(35-n-1));
	}
	PD.coord[na][0]+=f;
      }
      if (n==37) {
	  PD.coord[na][0]=pmflag*PD.coord[na][0];
      }
      if (n==38)
	pmflag=1.0;
      if (38<= n && n <42) {
	if (c==' ') {
	  f=0;
	}
	else if (c=='-') {
	  pmflag=-1.0;
	}
	else if (isdigit(c)) {
	  f=(double)(c-'0')*pow(10.0,(42-n-1));
	}
	PD.coord[na][1]+=f;
      }
      if (43<= n && n <46) {
	if (c==' ') {
	  f=0;
	}
	else if (isdigit(c)) {
	  f=(double)(c-'0')*pow(10.0,(43-n-1));
	}
	PD.coord[na][1]+=f;
      }
      if (n==46) {
	  PD.coord[na][1]=pmflag*PD.coord[na][1];
      }
      if (n==46)
	pmflag=1.0;
      if (46<= n && n <50) {
	if (c==' ') {
	  f=0;
	}
	else if (c=='-') {
	  pmflag=-1.0;
	}
	else if (isdigit(c)) {
	  f=(double)(c-'0')*pow(10.0,(50-n-1));
	}
	PD.coord[na][2]+=f;
      }
      if (51<= n && n <54) {
	if (c==' ') {
	  f=0;
	}
	else if (isdigit(c)) {
	  f=(double)(c-'0')*pow(10.0,(51-n-1));
	}
	PD.coord[na][2]+=f;
      }
      if (n==54) {
	  PD.coord[na][2]=pmflag*PD.coord[na][2];
      }
      PD.occupancy[na]=0.0;
      if (54<= n && n <57) {
	if (c==' ') {
	  f=0;
	}
	else if (isdigit(c)) {
	  f=(double)(c-'0')*pow(10.0,(57-n-1));
	}
	PD.occupancy[na]+=f;
      }
      if (58<= n && n <60) {
	if (c==' ') {
	  f=0;
	}
	else if (isdigit(c)) {
	  f=(double)(c-'0')*pow(10.0,(58-n-1));
	}
	PD.occupancy[na]+=f;
      }
      PD.tempfact[na]=0.0;
      if (60<= n && n <63) {
	if (c==' ') {
	  f=0;
	}
	else if (isdigit(c)) {
	  f=(double)(c-'0')*pow(10.0,(63-n-1));
	}
	PD.tempfact[na]+=f;
      }
      if (64<= n && n <=66) {
	if (c==' ') {
	  f=0;
	}
	else if (isdigit(c)) {
	  f=(double)(c-'0')*pow(10.0,(64-n-1));
	}
	PD.tempfact[na]+=f;
      }
    }
  }
  PD.na=na;
}

