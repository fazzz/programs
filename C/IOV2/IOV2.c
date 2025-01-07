
#include <stdio.h>

#include "IOV2.h"
#include "PT.h"

int IOV2_scantrj(FILE* trj,double ***crd,int MODE){
  int i=0,j=0,step=0;
  int c;
  int numatom;

  double f;

  numatom=AP.NATOM;

  while ((c=fscanf(trj,"%lf",&f))!=EOF){
    if ( MODE==AA ||
	 ( MODE==CA && strncmp(AP.IGRAPH[i],"CA",2)==0) ||
	 ( MODE==HV && strncmp(AP.IGRAPH[i],"H",1)!=0) ) {
      crd[i][j][step]=f;
    }
    ++j;
    if (j==3) {
      j=0;
      ++numatom;
    }
    if (i==numatom) {
      j=0;
      i=0;
      ++step;
      crd[i][j]=(double *)gcerealloc(crd[i][j],sizeof(double)*step);
    }
  }

  return step;
}

int IOV2_scantrj_wrange(FILE* trj,int inistep,int numstep,double ***crd,int MODE){
  int i=0,j=0,step=0;
  int c;
  int numatom;

  double f;

  numatom=AP.NATOM;

  while ((c=fscanf(trj,"%lf",&f))!=EOF){
    if (step > inistep) {
      if ( MODE==AA ||
	   ( MODE==CA && strncmp(AP.IGRAPH[i],"CA",2)==0) ||
	   ( MODE==HV && strncmp(AP.IGRAPH[i],"H",1)!=0) ) {
	crd[i][j][step]=f;
      }
    }
    ++j;
    if (j==3) {
      j=0;
      ++numatom;
    }
    if (i==numatom) {
      j=0;
      i=0;
      ++step;
      if (step > numstep) break;
      crd[i][j]=(double *)gcerealloc(crd[i][j],sizeof(double)*step);
    }
  }

  return step;
}
