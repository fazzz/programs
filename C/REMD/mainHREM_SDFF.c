
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PT.h"
#include "EF.h"
#include "REMD.h"

int main(int argc, char *argv[]) {
  int i,j;
  int numexchange,numreplica;
  int num,numatom;
  int numnb,num14;
  double temp;
  int waittime,numstep;

  char **crd,**top,*pn,*dirbase;
  char *minin,*rexin,*samin,*indexfilename;
  char crdname[40],topname[40];

  pn=getenv("pn");
  numreplica=atoi(getenv("numreplica"));
  numexchange=atoi(getenv("numexchange"));
  numstep=atoi(getenv("numstep"));
  numatom=atoi(getenv("numatom"));
  temp=atof(getenv("temp"));
  waittime=atoi(getenv("waittime"));
  dirbase=getenv("dirbase");
  numnb=atoi(getenv("numnb"));
  num14=atoi(getenv("num14"));
  minin=getenv("minin");
  rexin=getenv("rexin");
  samin=getenv("samin");
  indexfilename=getenv("indexfilename");

  crd=(char **)gcemalloc(sizeof(char *)*numreplica);
  top=(char **)gcemalloc(sizeof(char *)*numreplica);

  for (i=0;i<numreplica;++i) {
    sprintf(crdname,"crd%d",i+1);
    sprintf(topname,"top%d",i+1);
    crd[i]=getenv(crdname);
    top[i]=getenv(topname);
  }
  
  HREMD_vac_SDFF(numexchange,numreplica,minin,rexin,samin,crd,top,pn,dirbase,indexfilename,numnb,num14,temp,waittime,numstep,numatom);

  return 0;
}


