
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "CB.h"
#include "PT.h"
#include "TOPO.h"

int CB(double *crd,int MODE,double *bl) {
  int i,j,k;
  int numnb=0;
  double atom[2][3];

    if (MODE==INCH) {
    for (i=0;i<AP.NBONH;++i) {
      for (j=0;j<2;++j) 
	for (k=0;k<3;++k) atom[j][k]=crd[(AP.BH[i][j])+k];

      //      bl=(double *)gcerealloc(bl,sizeof(double)*(numnb+1));
      bl[numnb]=len(atom[0],atom[1]);
      ++numnb;
    }
  }
  for (i=0;i<AP.MBONA;++i) {
    for (j=0;j<2;++j) 
      for (k=0;k<3;++k) atom[j][k]=crd[(AP.BA[i][j])+k];
    
    //    bl=(double *)gcerealloc(bl,sizeof(double)*(numnb+1));
    bl[numnb]=len(atom[0],atom[1]);
    ++numnb;
  }

  return numnb;
}

