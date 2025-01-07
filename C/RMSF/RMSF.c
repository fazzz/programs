
#include "math.h"
#include "RMSF.h"
#include "mymath.h"

double *RMSF_prot(double ***crd/*atom*3*steps*/,int numatom,int numstep) {
  int i,j;
  double var[3];
  double *rmsf;
  
  rmsf=(double *)gcemalloc(sizeof(double)*numatom);

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) var[j]=MYMATH_var(crd[i][j],numstep);
    rmsf[i]=sqrt(var[0]*var[0]+var[1]*var[1]+var[2]*var[2]);
  }
  
  return rmsf;
}

