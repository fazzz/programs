
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PT.h"
#include "EF.h"
#include "REMD.h"

int main(int argc, char *argv[]) {
  int i,j;
  int numjudge,numreplica;
  int num;
  double temp;
  int waittime;

  char *dummy,*dummy2,*dummy3,*dummy4;

  char *crd,* vel,*top,*clt,*pn,*dirbase;
  char *inpmdini[20],*inpmdrst[20],*pepca[20];
  char inpininame[40],inprstname[40],pepcaname[40];

  dummy=getenv("numreplica");
  dummy2=getenv("numjudge");
  dummy3=getenv("tp");
  dummy4=getenv("waittime");

  numreplica=0;
  for (i=0;dummy[i]!=NULL;++i){
    num=dummy[i];numreplica=numreplica*10+num-'0';
  }
  numjudge=0;
  for (i=0;dummy2[i]!=NULL;++i){
    num=dummy2[i];numjudge=numjudge*10+num-'0';
  }
  temp=0.0;
  for (i=0;dummy3[i]!=NULL;++i){
    num=dummy3[i];temp=temp*10+(double)(num-'0');
  }
  waittime=0;
  for (i=0;dummy4[i]!=NULL;++i){
    num=dummy4[i];waittime=waittime*10+num-'0';
  }

  crd=getenv("crd");
  vel=getenv("vel");
  top=getenv("top");
  clt=getenv("clt");
  pn=getenv("pn");
  dirbase=getenv("dirbase");
  for (i=0;i<numreplica;++i) {
    sprintf(inpininame,"inpmdini%d",i+1);
    sprintf(inprstname,"inpmdrst%d",i+1);
    sprintf(pepcaname, "pepca%d"   ,i+1);
    inpmdini[i]=getenv(inpininame);
    inpmdrst[i]=getenv(inprstname);
    pepca[i]=getenv(pepcaname);
  }
  
  HreplicaExchange(numjudge,numreplica,crd,vel,top,clt,pn,inpmdini,inpmdrst,pepca,dirbase,temp,waittime);
  
  return 0;
}


