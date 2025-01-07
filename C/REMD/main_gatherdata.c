
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
  int numstep;
  char *pn,*dirbase;

  pn=getenv("pn");
  numreplica=atoi(getenv("numreplica"));
  numexchange=atoi(getenv("numexchange"));
  numstep=atoi(getenv("numstep"));
  numatom=atoi(getenv("numatom"));
  dirbase=getenv("dirbase");

  gatherdata(numexchange,numreplica,numstep,numatom,dirbase,pn);

  return 0;
}


