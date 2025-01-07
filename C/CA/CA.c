
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "CA.h"
#include "PT.h"
#include "TOPO.h"

int CA(double *crd,int MODE,double *ba) {
  int i,j,k;
  int numna=0;
  double atom[3][3];

    if (MODE==INCH) {
      for (i=0;i<AP.NTHETH;++i) {
	for (j=0;j<2;++j) 
	  for (k=0;k<3;++k) atom[j][k]=crd[(AP.TH[i][j])+k];

	//	ba=(double *)gcerealloc(ba,sizeof(double)*(numna+1));
	ba[numna]=ang(atom[0],atom[1],atom[2]);
	++numna;
      }
    }
  for (i=0;i<AP.MTHETA;++i) {
    for (j=0;j<3;++j) 
      for (k=0;k<3;++k) atom[j][k]=crd[(AP.TA[i][j])+k];
    
    //    ba=(double *)gcerealloc(ba,sizeof(double)*(numna+1));
    ba[numna]=ang(atom[0],atom[1],atom[2]);
    ++numna;
  }

  return numna;
}

