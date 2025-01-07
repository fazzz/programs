
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PT.h"
#include "EF.h"
#include "REMD.h"

int main(int argc, char *argv[]) {
  int i,j,numflag;
  int numjudge,numreplica;
  int num;
  
  char *dummy,*dummy2;

  char *crd,* vel,*top,*clt,*pn,*dirbase;
  char *inpmdini[20],*inpmdrst[20],*pepca[20];
  char inpininame[40],inprstname[40],pepcaname[40];

  dummy=getenv("numreplica");
  dummy2=getenv("numjudge");

  numreplica=0;
  for (i=0;dummy[i]!=NULL;++i){
    num=dummy[i];numreplica=numreplica*10+num-'0';
  }
  numjudge=0;
  for (i=0;dummy2[i]!=NULL;++i){
    num=dummy2[i];numjudge=numjudge*10+num-'0';
  }

  crd=getenv("crd");
  vel=getenv("vel");
  top=getenv("top");
  clt=getenv("clt");
  pn=getenv("pn");
  dirbase=getenv("dirbase");
  for (i=0;i<10;++i) {
    sprintf(inpininame,"inpmdini%d",i+1);
    sprintf(inprstname,"inpmdrst%d",i+1);
    sprintf(pepcaname, "pepca%d"   ,i+1);
    inpmdini[i]=getenv(inpininame);
    inpmdrst[i]=getenv(inprstname);
    pepca[i]=getenv(pepcaname);
  }
  
  TreplicaExchange(numjudge,numreplica,crd,vel,top,clt,pn,inpmdini,inpmdrst,pepca,dirbase);
  
  return 0;
}


